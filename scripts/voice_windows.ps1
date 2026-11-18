param([Parameter(Mandatory=$true)][string]$OutputDirectory)
$ErrorActionPreference = 'Stop'
Add-Type -AssemblyName System.Speech
$lines = @(
'The next step is three-dimensional immersion cooling.',
'A heated chip below. A cold plate above.',
'Built on a Cartesian grid.',
'We reconstruct the moving liquid and vapor interface in three dimensions.',
'Temperature and velocity come directly from the calculation.',
'Early numerical evolution. The full cycle remains under study.'
)
$cuts = @(0,3,7,10,15,20,24)
New-Item -ItemType Directory -Path $OutputDirectory -ErrorAction Stop | Out-Null
$speaker = New-Object System.Speech.Synthesis.SpeechSynthesizer
try {
  $available = @($speaker.GetInstalledVoices() | Where-Object {$_.Enabled -and $_.VoiceInfo.Culture.Name -like 'en-*'})
  if ($available.Count -eq 0) { throw 'No local English synthesis voice installed.' }
  $speaker.SelectVoice($available[0].VoiceInfo.Name)
  $records = @()
  for ($i=0; $i -lt $lines.Count; $i++) {
    $path = Join-Path $OutputDirectory ('VOICE_{0:D2}.wav' -f $i)
    $speaker.Rate = 4
    $speaker.SetOutputToWaveFile($path)
    $speaker.Speak($lines[$i])
    $speaker.SetOutputToNull()
    $records += [pscustomobject]@{scene=$i; start=$cuts[$i]; end=$cuts[$i+1]; text=$lines[$i]; nativeRate=4; postTimeStretch=$false}
  }
  [pscustomobject]@{voice=$available[0].VoiceInfo.Name; synthetic=$true; scenes=$records} | ConvertTo-Json -Depth 5 | Set-Content -LiteralPath (Join-Path $OutputDirectory 'VOICE.json') -Encoding utf8
} finally { $speaker.Dispose() }
