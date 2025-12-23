Param(
  [string]$Exe = ".\b02_executable.exe",
  [string]$Geometry = "primitive",
  [string]$Out = "outputs",
  [string]$Tag = "",
  [int]$Events = 1000000,
  [int]$Chunks = 1,
  [int]$Seed1 = 12345,
  [int]$Seed2 = 67890,
  [switch]$MakePaperOutputs
)

if (-not $Tag) {
  $Tag = "run_" + (Get-Date -Format "yyyyMMdd_HHmmss")
}

$OutDir = Join-Path $Out $Tag
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
$EventsPerChunk = [int][Math]::Ceiling($Events / [double]$Chunks)
$chunkFiles = @()

for ($i=0; $i -lt $Chunks; $i++) {
  $chunkTag = "{0}_chunk{1:d2}" -f $Tag, $i
  $outfile = Join-Path $OutDir ("B02ntuples_{0}.root" -f $chunkTag)
  $logfile = Join-Path $OutDir ("run_{0}.log" -f $chunkTag)
  $macro = Join-Path $OutDir ("macro_{0}.mac" -f $chunkTag)
  $seed1 = $Seed1 + $i
  $seed2 = $Seed2 + $i
  @"
/analysis/setFileName $outfile
/analysis/ntuple/setFileName 0 $outfile
/analysis/ntuple/setFileName 1 $outfile
/runAction/useTimeSeed false
/runAction/useFixedSeeds true
/runAction/seed1 $seed1
/runAction/seed2 $seed2
/run/printProgress 10000
/generator/SmithActivation true
/generator/useFixedEnergy false
/run/initialize
/run/beamOn $EventsPerChunk
"@ | Set-Content -Path $macro -Encoding ASCII

  Write-Host "[run_simccd_1M] chunk $i : launching $Exe"
  & $Exe --geometry $Geometry --no-vis $macro *> $logfile
  if ($LASTEXITCODE -ne 0) {
    Write-Host "[run_simccd_1M] chunk $i failed; see $logfile"
    exit $LASTEXITCODE
  }
  if (-not (Test-Path $outfile)) {
    Write-Host "[run_simccd_1M] chunk $i produced no ROOT file ($outfile)"
    exit 1
  }
  $chunkFiles += $outfile
}

$merged = ""
if ($chunkFiles.Count -gt 1) {
  if (Get-Command hadd -ErrorAction SilentlyContinue) {
    $merged = Join-Path $OutDir ("B02ntuples_{0}_merged.root" -f $Tag)
    Write-Host "[run_simccd_1M] merging chunks into $merged"
    hadd -f $merged @chunkFiles | Out-Null
  } else {
    Write-Host "[run_simccd_1M] hadd not found; skipping merge."
  }
} else {
  $merged = $chunkFiles[0]
}

@{
  exe = $Exe
  geometry = $Geometry
  events = $Events
  chunks = $Chunks
  events_per_chunk = $EventsPerChunk
  seed1 = $Seed1
  seed2 = $Seed2
  output_dir = $OutDir
  merged_file = $merged
} | ConvertTo-Json | Set-Content -Path (Join-Path $OutDir "run_metadata.json") -Encoding ASCII

if ($MakePaperOutputs -and $merged) {
  python analysis/make_paper_outputs.py --input $merged --output paper_outputs --tag $Tag
}

Write-Host "[run_simccd_1M] DONE. Outputs in $OutDir"
