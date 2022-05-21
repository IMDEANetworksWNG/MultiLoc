# Empty the data file
echo "" > /tmp/aoa_measurements.txt

# Reset the interface
ip link set dev wlan0 down
ip link set dev wlan0 up

# Reset the supplicant
killall wpa_supplicant

# Connect 
wpa_supplicant -D nl80211  -i wlan0 -c /etc/wpa_supplicant.conf -B

# This will help to check if we cannot connect
i=0
connected=0

while true ; do  # Loop until connected

    connected=$(cat /sys/kernel/debug/ieee80211/phy0/wil6210/stations | grep connected | wc -l)
                                             
    if [[ $connected -eq 1 ]]; then
      break
    fi

    i=$((i+1))

    # Break after a certain number of measurements
    if [[ $i -eq 10 ]]; then
      echo "[AoA] We cannot connect after 10 seconds"
      echo "Unreachable" > /tmp/aoa_measurements.txt
      exit
    fi

    sleep 1
done

echo "[AoA] We are connected now"

# A counter
i=0

while true ; do  # Loop until interval has elapsed.

    connected=$(cat /sys/kernel/debug/ieee80211/phy0/wil6210/stations | grep connected | wc -l)

    # Get AoA
    echo -n -e '\x04\xd6\xaa\x6e\xcc\x9f' | iw dev wlan0 vendor recv 0x001374 0x93 -

    # Save the measurement
    echo $(dmesg | tail -n1) >> /tmp/aoa_measurements.txt

    i=$((i+1))

    # Break when no connected
    if [[ $connected -eq 0 ]]; then
      echo "[AoA] We got " $i "measurements" 
      break
    fi

    # Break after a certain number of measurements
    if [[ $i -eq 300 ]]; then
      echo "[AoA] We got 300 measurements"
      break
    fi

done
