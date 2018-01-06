Do {
    If (!(Get-Process -name Client -ErrorAction SilentlyContinue)) {
           Write-Host "Wating for Process to Start"
           Start-Sleep -Seconds 2.0
		   Start-Process "C:\Projects\Analytics\AnomalyDetection\Client\Data\client.exe"
    }Else {
           Write-Host 'Process has Started'
           While (Get-Process -name Client -ErrorAction SilentlyContinue) {
                    Write-Host 'Waiting for Process to be Stopped'
                    Start-Sleep -Seconds 2.0
           }
           Write-Host 'Process Stopped' #; $Status = 'Done'
    }
}Until ($Status)s