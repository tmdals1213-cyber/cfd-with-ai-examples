from pathlib import Path
import json,io,subprocess,argparse,hashlib
import numpy as np
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image,ImageDraw,ImageFont
ROOT=Path(__file__).resolve().parents[1]
BG='#08131f';FG='#f0f5f8';TEAL='#4dddc6';ORANGE='#f2a352';GRAY='#91a9ba'
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--ffmpeg',default='ffmpeg');a=ap.parse_args();entry=json.loads((ROOT/'examples/legacy_catalog.json').read_text())['15_runner_final'];case=ROOT/entry['directory'];d=json.loads((case/'case.json').read_text());q=json.loads((case/'RECORDED_PRESSURE_SUMMARY.json').read_text());fontdir=Path(matplotlib.get_data_path())/'fonts/ttf'
 def font(n,b=False):return ImageFont.truetype(str(fontdir/('DejaVuSans-Bold.ttf'if b else'DejaVuSans.ttf')),n)
 fig,ax=plt.subplots(figsize=(8,8),facecolor=BG);ax.set_facecolor(BG);positions=np.array(d['offsets_m'])
 # Plan view: forward +X is up; lateral Z is horizontal. Every marker is an actual prescribed offset.
 for i,(x,y,z)in enumerate(positions):
  c=ORANGE if i==0 else TEAL
  ax.scatter(z,x,s=650,color=c,edgecolor=FG,linewidth=.6)
 ax.annotate('AIR',xy=(-2.3,.3),xytext=(-2.3,3.4),color=GRAY,fontsize=17,ha='center',arrowprops={'arrowstyle':'->','color':GRAY,'lw':2})
 ax.text(0,-1.35,'Protected runner',ha='center',color=ORANGE,fontsize=20);ax.set(xlim=(-3,3),ylim=(-1.65,4.25));ax.set_aspect('equal');ax.axis('off');fig.tight_layout(pad=.3);b=io.BytesIO();fig.savefig(b,format='png',dpi=160,facecolor=BG);plt.close(fig);b.seek(0);formation=Image.open(b).convert('RGB')
 fig,ax=plt.subplots(figsize=(9,6),facecolor=BG);ax.set_facecolor(BG);values=[q['meanSolo_N'],q['meanPaced_N']];bars=ax.barh([1,0],values,color=[ORANGE,TEAL],height=.45);ax.axvline(0,color=GRAY,lw=1);ax.set_yticks([1,0],['SOLO','PACED'],color=FG,fontsize=19);ax.set_xlim(-1.5,10);ax.set_xticks([0,4,8],['0','4','8'],color=GRAY,fontsize=16);ax.set_xlabel('Pressure component [N]',color=FG,fontsize=19);ax.spines[['top','right','left']].set_visible(False);ax.spines['bottom'].set_color(GRAY)
 for y,v in zip([1,0],values):ax.text(max(v,0)+.25,y,f'{v:+.3f} N',va='center',color=FG,fontsize=20)
 fig.tight_layout();b=io.BytesIO();fig.savefig(b,format='png',dpi=150,facecolor=BG);plt.close(fig);b.seek(0);bar=Image.open(b).convert('RGB')
 records=[]
 for portrait in [False,True]:
  W,H=(1080,1920)if portrait else(1920,1080);form='PORTRAIT'if portrait else'LANDSCAPE';slides=[]
  titles=['RUNNER DRAFTING','SAME SPEED. SAME MOTION.','PRESSURE COMPONENT','A CONTROLLED COMPARISON'];subtitles=['Solo vs seven pacemakers','21.1 km/h  /  192 spm','Recorded mean: 1.25–2.50 s','Pressure-only diagnostic']
  for i in range(4):
   im=Image.new('RGB',(W,H),BG);draw=ImageDraw.Draw(im);title=titles[i]
   if portrait:title=title.replace(' DRAFTING','\nDRAFTING').replace('SAME SPEED. SAME MOTION.','SAME SPEED.\nSAME MOTION.').replace(' COMPONENT','\nCOMPONENT').replace('A CONTROLLED COMPARISON','A CONTROLLED\nCOMPARISON')
   draw.text((W*.07,H*.04),'RUNNER AERODYNAMICS',font=font(24,True),fill=TEAL);draw.multiline_text((W*.07,H*.12),title,font=font(61 if portrait else 64,True),fill=FG,spacing=10)
   asset=formation if i<2 else bar;asset=asset.copy();rect=(W*.05,H*.31,W*.9,H*.43)if portrait else(W*.40,H*.25,W*.56,H*.58);x,y,w,h=map(int,rect);asset.thumbnail((w,h),Image.Resampling.LANCZOS);im.paste(asset,(x+(w-asset.width)//2,y+(h-asset.height)//2))
   if portrait:draw.text((W*.07,H*.80),subtitles[i],font=font(37,True),fill=ORANGE)
   else:draw.text((W*.07,H*.50),subtitles[i],font=font(32,True),fill=ORANGE)
   footer='Recorded formation and means'if i<2 else'Viscous stress excluded; exterior traction unqualified'
   draw.text((W*.07,H*.93),footer,font=font(24 if portrait else 23),fill=GRAY);slides.append(im)
  slides[2].save(ROOT/'media'/('RUNNER_POSTER_'+form+'.png'));target=ROOT/'media'/('RUNNER_RECORDED_SUMMARY_'+form+'.mp4');assert not target.exists()
  command=[a.ffmpeg,'-hide_banner','-loglevel','error','-n','-f','rawvideo','-pix_fmt','rgb24','-s',f'{W}x{H}','-r','24','-i','pipe:0','-an','-c:v','libx264','-threads','2','-crf','20','-preset','veryfast','-pix_fmt','yuv420p','-movflags','+faststart',str(target)];p=subprocess.Popen(command,stdin=subprocess.PIPE)
  try:
   for slide in slides:
    data=slide.tobytes()
    for j in range(96):p.stdin.write(data)
   p.stdin.close();assert p.wait(timeout=60)==0
  finally:
   if p.poll()is None:p.kill();p.wait()
  subprocess.run([a.ffmpeg,'-hide_banner','-loglevel','error','-i',str(target),'-f','null','-'],check=True,timeout=60)
  records.append({'file':'media/'+target.name,'seconds':16,'fps':24,'width':W,'height':H,'sha256':hashlib.sha256(target.read_bytes()).hexdigest(),'kind':'New presentation of recorded positions and historical pressure-component means; not recovered flow footage','newCFDStatesGenerated':False,'physicalTotalDragQualified':False,'fullDecodePassed':True});print(json.dumps(records[-1]),flush=True)
 (ROOT/'provenance/runner_summary_media.json').write_text(json.dumps(records,indent=2)+'\n')
if __name__=='__main__':main()
