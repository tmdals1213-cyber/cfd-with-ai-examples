from pathlib import Path
import argparse,json,hashlib,subprocess,io
from PIL import Image
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser();p.add_argument('--ffmpeg',default='ffmpeg');a=p.parse_args();parts=[('VEHICLE_CAR10_LANDSCAPE.mp4',1),('RUNNER_RECORDED_SUMMARY_LANDSCAPE.mp4',4),('IMMERSION_COOLING_LANDSCAPE_DRAFT01.mp4',10)];target=ROOT/'media/PROJECTS_OVERVIEW.mp4';assert not target.exists();args=[];filters=[];records=[]
 for i,(name,start)in enumerate(parts):
  source=ROOT/'media'/name;args+=['-i',str(source)];filters.append(f'[{i}:v]trim=start={start}:duration=8,setpts=PTS-STARTPTS,scale=1920:1080,setsar=1,fps=24[v{i}]');records.append({'source':name,'sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'presentationStart_s':start,'presentationDuration_s':8})
 filters.append('[v0][v1][v2]concat=n=3:v=1:a=0[out]');subprocess.run([a.ffmpeg,'-hide_banner','-loglevel','error','-n',*args,'-filter_complex',';'.join(filters),'-map','[out]','-an','-c:v','libx264','-threads','2','-crf','23','-preset','fast','-pix_fmt','yuv420p','-movflags','+faststart',str(target)],check=True,timeout=180)
 subprocess.run([a.ffmpeg,'-hide_banner','-loglevel','error','-i',str(target),'-f','null','-'],check=True,timeout=90)
 sheet=Image.new('RGB',(960,540*3),'#08131f')
 for i,t in enumerate([4,12,20]):
  b=subprocess.check_output([a.ffmpeg,'-hide_banner','-loglevel','error','-ss',str(t),'-i',str(target),'-frames:v','1','-vf','scale=960:540','-f','image2pipe','-c:v','png','pipe:1']);im=Image.open(io.BytesIO(b)).convert('RGB');sheet.paste(im,(0,540*i))
 sheet.save(ROOT/'media/PROJECTS_CONTACT.png');d={'file':'media/'+target.name,'seconds':24,'fps':24,'components':records,'newCFDStatesGenerated':False,'sourceKinds':['Single final-field car camera orbit','Recorded runner offsets and pressure means','Actual1–20ms XYZ PLIC/field showcase'],'fullDecodePassed':True,'sha256':hashlib.sha256(target.read_bytes()).hexdigest()};(ROOT/'provenance/projects_overview.json').write_text(json.dumps(d,indent=2)+'\n');print(json.dumps({'file':target.name,'bytes':target.stat().st_size,'fullDecodePassed':True}))
if __name__=='__main__':main()
