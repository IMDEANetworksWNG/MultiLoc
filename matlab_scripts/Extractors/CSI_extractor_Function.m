function [cmplxall_raw_all] = CSI_extractor_Function(FILE)

%  = strcat("pcap_files/",index_set(id_set),"/./trace_", string(id_meas),".pcap");

HOFFSET = 16;           % header offset
BW = 80;
NFFT = BW*3.2;          % fft size
p = readpcap();
p.open(FILE);
n = length(p.all());
p.from_start();
k = 1;

prevfwcnt = -1;
nss_config = 4;
rxcore_config = 4;
mask_toprocess = nss_config * rxcore_config;
slice = uint32(0);

processedmask = mask_toprocess;
processedpacket = 0;
chars='|/-\';
output = '';
fprintf(1, '.');

prevfwmask = -1;

csi_store = {};
packets = 1;

while (k <= n)
for jj = 1:length(output) + 1, fprintf(1, '\b'); end;
charjj = mod(k, length(chars)) + 1;
output = sprintf(' %.0f/100', k / n * 100);
fprintf(1, '%c%s', chars(charjj), output);
f = p.next();
if isempty(f)
    disp('no more frames');
    return;
end
if f.header.orig_len-(HOFFSET-1)*4 < NFFT*4
    disp('skipped frame with incorrect size');
    return;
end

% extract number of packet processed by the firmware and rxcore for this packet
fwcnt = bitand(f.payload(15), 65535 * 65536) / 65536;
fwmask = bitand(f.payload(14), 255 * 65536) / 65536;
rxcore = double(bitand(fwmask, 12) / 4 + 1);
nss = double(bitand(fwmask, 3) + 1);
% disp(sprintf('%d %d %d %d', [fwcnt fwmask rxcore nss]));

rxcore = double(rxcore);

% if fwcnt wrapped around, then this is like when we have a new packet
realprevfwcnt = prevfwcnt;
if fwcnt < prevfwcnt,
  disp 'fwcnt wrapped around...';
  prevfwcnt = fwcnt - 1;
end;

% if this data is for a new packet, reset number of rxcore for this packet
% report if for the previous packet we did not received all rxcore data
if fwcnt > prevfwcnt,
  if processedmask < mask_toprocess,
    disp(sprintf('\nmissing data for packet %d\n', realprevfwcnt));
    for jj = 1:length(output), fwrite(1, 10); end;
    % drawnow('update');
  end;
  processedmask = 0;
  prevfwcnt = fwcnt;
  prevfwmask = -1;
  slice = uint32(0);
end;

processedmask = processedmask + 1;

% disp(sprintf('%d %d %d %d %d', [fwmask, processedmask, mask_toprocess, prevfwcnt, fwcnt]));

% this should not happen
if processedmask > mask_toprocess,
  disp(sprintf('More than %d masks for this packet, terminating...', mask_toprocess));
  break;
end;

if size(slice, 1) == 1,
  slice = uint32(zeros([length(f.payload) 4]));
end;

slice(:, fwmask + 1) = f.payload;

if fwmask < prevfwmask,
  disp(sprintf('\nfwmask goes backward, skipping CSI %d\n', fwcnt));
  for jj = 1:length(output), fwrite(1, 10); end;
  % drawnow('update');
  processedmask = 0;
  prevfwcnt = fwcnt;
  prevfwmask = -1;
  slice = uint32(0);
end;

prevfwmask = fwmask;

% if we have enough data for this packet process everything
if processedmask == mask_toprocess,
  cmplxall = zeros([256 mask_toprocess]);
  cmplxall_raw = zeros([256 mask_toprocess]);

  % extract CSI
  % 1) extraction
  for jj = 1:mask_toprocess,
    payload = slice(:, jj);
    H = payload(HOFFSET:HOFFSET+NFFT-1);
    Hout = unpack_float(int32(1), int32(NFFT), H);
    Hout = reshape(Hout,2,[]).';
    cmplx = double(Hout(1:NFFT,1))+1j*double(Hout(1:NFFT,2));
    cmplx = cmplx([129:256 1:128]);

    cmplxall_raw(:, jj) = cmplx;
% normalisation
    % cmplx = cmplx ./ reference(:, jj);
    cmplxall(:, jj) = cmplx;

  end;

  cmplxall_raw_all(packets,:,:) = cmplxall_raw;


  packets = packets + 1;

  if 0,
    figure(1); clf;
    for jj = 1:mask_toprocess,
      subplot(nss_config, rxcore_config, jj);
      plot(abs(cmplxall_raw(:, jj)));
    end;
    drawnow;
  end;

end;

k = k + 1;
end


packets = packets - 1;

end
