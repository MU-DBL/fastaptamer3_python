#!/bin/bash
# setup-docker-resources.sh
# Configures Docker Desktop on macOS to use more CPU/memory for large SELEX file processing.
# Run once before first use. No admin rights required.

SETTINGS="$HOME/Library/Group Containers/group.com.docker/settings.json"

# Detect system resources
TOTAL_RAM_MB=$(( $(sysctl -n hw.memsize) / 1024 / 1024 ))
TOTAL_CPUS=$(sysctl -n hw.logicalcpu)

# Allocate 75%, with sane minimums
DOCKER_RAM_MB=$(( TOTAL_RAM_MB * 75 / 100 ))
DOCKER_CPUS=$(( TOTAL_CPUS * 75 / 100 ))
[ $DOCKER_RAM_MB -lt 4096 ] && DOCKER_RAM_MB=4096
[ $DOCKER_CPUS -lt 2 ] && DOCKER_CPUS=2

DOCKER_RAM_GB=$(( DOCKER_RAM_MB / 1024 ))

echo ""
echo "System detected:"
echo "  Total RAM : $(( TOTAL_RAM_MB / 1024 )) GB"
echo "  Total CPUs: $TOTAL_CPUS"
echo ""
echo "Docker Desktop will be configured to use:"
echo "  Memory : ${DOCKER_RAM_GB} GB  (75% of total)"
echo "  CPUs   : $DOCKER_CPUS  (75% of total)"
echo ""

if [ ! -f "$SETTINGS" ]; then
    echo "Error: Docker Desktop settings not found at:"
    echo "  $SETTINGS"
    echo "Make sure Docker Desktop is installed and has been run at least once."
    exit 1
fi

# Backup existing settings
cp "$SETTINGS" "${SETTINGS}.bak"
echo "Existing settings backed up to ${SETTINGS}.bak"

# Update settings using Python (available on all Macs)
python3 - <<EOF
import json

path = "$SETTINGS"
with open(path) as f:
    settings = json.load(f)

settings['cpus'] = $DOCKER_CPUS
settings['memoryMiB'] = $DOCKER_RAM_MB

with open(path, 'w') as f:
    json.dump(settings, f, indent=2)

print("Settings updated successfully.")
EOF

echo ""
echo "Please restart Docker Desktop to apply changes:"
echo "  Docker Desktop menu → Restart"
echo ""
