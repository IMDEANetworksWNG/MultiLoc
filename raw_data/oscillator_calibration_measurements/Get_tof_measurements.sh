# Empty the file
echo "" > /tmp/tof_measurements.txt

# Total number of measurements
total_measurements=0
connected=0

# Clear dmesg
dmesg -c > /dev/null 2>&1

while true; do

    killall wpa_supplicant

    ip link set dev wlan0 down; ip link set dev wlan0 up

    wpa_supplicant -D nl80211  -i wlan0 -c /etc/wpa_supplicant.conf -B

    retries=0

    while true ; do  # Loop until connected

        connected=$(cat /sys/kernel/debug/ieee80211/phy0/wil6210/stations | grep connected | wc -l)
                                                                       
        if [[ $connected -eq 1 ]]; then
          break
        fi

        retries=$((retries+1))
        sleep 1

        if [[ $retries -eq 10 ]]; then

            break
        fi
    done

    # Check again to escape from outer loop
    if [[ $retries -eq 10 ]]; then

        echo "There is no connection"
        echo "Unreachable" >> /tmp/tof_measurements.txt
        break
    fi

    if [[ $connected -eq 1 ]]; then
        echo "We are connected now"

        # ToF command
        echo -n -e '\x04\xd6\xaa\x6e\xcc\x9f' | iw dev wlan0 vendor recv 0x001374 0x81 -

        # Check if there are any measurements
        num_measurements=$(dmesg | grep Measurement | wc -l )

        if [[ $num_measurements -gt 0 ]]; then
      
            # Save the measurements
            echo "$(dmesg -c | grep Measurement)" >> /tmp/tof_measurements.txt

            # Increment the number of measurements
            total_measurements=$((total_measurements+num_measurements))

            echo $num_measurements" new measurements ("$total_measurements"/100)" 
        fi

        # Break after a certain number of measurements
        if [[ $total_measurements -ge 100 ]]; then
          echo "We got "$total_measurements" measurements total"
          break
        fi
    fi
done

killall wpa_supplicant
