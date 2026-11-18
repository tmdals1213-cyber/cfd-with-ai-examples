from pathlib import Path
import argparse,subprocess,json,hashlib
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 p=argparse.ArgumentParser();p.add_argument('--outputs',type=Path,required=True);p.add_argument('--ffmpeg',default='ffmpeg');a=p.parse_args()
 selections=[('car10_english_overlay_film_20261102/final/CAR10_English_2560x1440.mp4','VEHICLE_CAR10_LANDSCAPE.mp4',1920,1080,'Camera orbit around the same final field; not unsteady time evolution.'),('car10_english_overlay_film_20261102/final/CAR10_English_SNS_1080x1920.mp4','VEHICLE_CAR10_PORTRAIT.mp4',1080,1920,'Same-time final-field comparison with camera motion.'),('drivaer_unsteady_demo/car_wake_preview.mp4','VEHICLE_DRIVAER_WAKE.mp4',1280,720,'Historical actual saved wake preview; model/grid limits retained.'),('ev_aero_comparison/drag_comparison/final/Model3_vs_IONIQ5_drag_9x16_final.mp4','VEHICLE_EV_PAIR.mp4',720,1280,'Historical qualitative artist-model EV comparison.'),('suv_aero_comparison/drag_comparison/final/URUS_vs_GCLASS_drag_9x16_final.mp4','VEHICLE_SUV_PAIR.mp4',720,1280,'Historical qualitative artist-model SUV comparison.')]
 records=[]
 for relative,name,w,h,scope in selections:
  source=a.outputs/relative;assert source.is_file();target=ROOT/'media'/name;assert not target.exists()
  cmd=[a.ffmpeg,'-hide_banner','-loglevel','error','-n','-i',str(source),'-vf',f'scale={w}:{h}:force_original_aspect_ratio=decrease,pad={w}:{h}:(ow-iw)/2:(oh-ih)/2','-r','24','-c:v','libx264','-threads','2','-preset','fast','-crf','26','-pix_fmt','yuv420p','-c:a','aac','-b:a','128k','-movflags','+faststart',str(target)]
  subprocess.run(cmd,check=True,timeout=180);subprocess.run([a.ffmpeg,'-hide_banner','-loglevel','error','-i',str(target),'-f','null','-'],check=True,timeout=90)
  records.append({'file':'media/'+name,'originalName':relative,'originalSHA256':sha(source),'sha256':sha(target),'bytes':target.stat().st_size,'scope':scope,'newCFDStatesGenerated':False,'fullDecodePassed':True})
  print(json.dumps(records[-1]),flush=True)
 (ROOT/'provenance/vehicle_media.json').write_text(json.dumps(records,indent=2)+'\n')
if __name__=='__main__':main()
