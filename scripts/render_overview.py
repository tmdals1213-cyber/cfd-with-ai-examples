from pathlib import Path
import argparse,subprocess,json,hashlib
from PIL import Image,ImageDraw,ImageFont
import matplotlib
ROOT=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser();p.add_argument('--ffmpeg',default='ffmpeg');p.add_argument('--output',type=Path,required=True);a=p.parse_args()
 slides=[('NATURAL CONVECTION','10 W · 10 s · EOS + CHT','water_cht_10s.png'),('A PHASE-CHANGE REFERENCE','Planar Stefan · one grid','stefan_reference.png'),('HEAT-DRIVEN VAPOR','100 W · 0.480 s · growth and rise','heated_bubble_480ms.png'),('A COLD CEILING','2D dielectric pool · early transient','dielectric_2d_time_comparison.png'),('INTO THREE DIMENSIONS','512,000 cells · 20 ms · 8 MPI','immersion_3d_20ms.png')]
 fontdir=Path(matplotlib.get_data_path())/'fonts/ttf'
 def font(n):return ImageFont.truetype(str(fontdir/'DejaVuSans.ttf'),n)
 frames=[]
 for title,sub,name in slides:
  im=Image.new('RGB',(1920,1080),'#08131f');d=ImageDraw.Draw(im);d.text((90,55),title,font=font(54),fill='white');d.text((90,130),sub,font=font(29),fill='#4dddc6');q=Image.open(ROOT/'media/previews'/name).convert('RGB');q.thumbnail((1740,770),Image.Resampling.LANCZOS);im.paste(q,((1920-q.width)//2,210+(770-q.height)//2));d.text((90,1020),'RECORDED EXAMPLES  /  NUMERICAL DEMONSTRATIONS',font=font(20),fill='#92a9b9');frames.append(im)
 cmd=[a.ffmpeg,'-hide_banner','-loglevel','error','-n','-f','rawvideo','-pix_fmt','rgb24','-s','1920x1080','-r','24','-i','pipe:0','-an','-c:v','libx264','-threads','2','-crf','20','-preset','veryfast','-pix_fmt','yuv420p','-movflags','+faststart',str(a.output)]
 proc=subprocess.Popen(cmd,stdin=subprocess.PIPE)
 try:
  for im in frames:
   raw=im.tobytes()
   for i in range(96):proc.stdin.write(raw)
  proc.stdin.close();assert proc.wait(timeout=60)==0
 finally:
  if proc.poll()is None:proc.kill();proc.wait()
 subprocess.run([a.ffmpeg,'-hide_banner','-loglevel','error','-i',str(a.output),'-f','null','-'],check=True,timeout=60)
 print(json.dumps({'file':a.output.name,'bytes':a.output.stat().st_size,'seconds':20,'fps':24,'fullDecodePassed':True,'kind':'Historical result slide summary, no new CFD states generated','sha256':hashlib.sha256(a.output.read_bytes()).hexdigest()}))
if __name__=='__main__':main()
