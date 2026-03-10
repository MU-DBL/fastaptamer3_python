# setup-docker-resources.ps1
# Configures WSL2 to give Docker more CPU/memory for large SELEX file processing.
# Run once as a regular user (no admin required).

$wslConfigPath = "$env:USERPROFILE\.wslconfig"

# Detect system resources
$totalRamGB  = [math]::Round((Get-CimInstance Win32_ComputerSystem).TotalPhysicalMemory / 1GB)
$totalCPUs   = (Get-CimInstance Win32_Processor | Measure-Object -Property NumberOfLogicalProcessors -Sum).Sum

# Allocate 75% of RAM and CPUs to WSL2, with sane minimums
$wslRamGB    = [math]::Max(4, [math]::Floor($totalRamGB * 0.75))
$wslCPUs     = [math]::Max(2, [math]::Floor($totalCPUs * 0.75))
$swapGB      = [math]::Max(4, [math]::Floor($wslRamGB / 2))

Write-Host ""
Write-Host "System detected:"
Write-Host "  Total RAM : ${totalRamGB} GB"
Write-Host "  Total CPUs: $totalCPUs"
Write-Host ""
Write-Host "WSL2 will be configured to use:"
Write-Host "  Memory : ${wslRamGB} GB  (75% of total)"
Write-Host "  CPUs   : $wslCPUs  (75% of total)"
Write-Host "  Swap   : ${swapGB} GB"
Write-Host ""

$config = @"
[wsl2]
memory=${wslRamGB}GB
processors=$wslCPUs
swap=${swapGB}GB
"@

# Backup existing config if present
if (Test-Path $wslConfigPath) {
    $backup = "$wslConfigPath.bak"
    Copy-Item $wslConfigPath $backup
    Write-Host "Existing .wslconfig backed up to $backup"
}

Set-Content -Path $wslConfigPath -Value $config
Write-Host ".wslconfig written to $wslConfigPath"
Write-Host ""
Write-Host "Restarting WSL2 to apply changes..."
wsl --shutdown
Write-Host "Done. Please restart Docker Desktop if it is running."
Write-Host ""
