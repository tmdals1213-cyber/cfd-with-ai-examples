"""Black logo/drag-chart presentation; every model stays visible in overview."""
import json,math,sys
from pathlib import Path
H=Path(__file__).resolve().parent;R=H.parents[1];O=R/'outputs/car10_final_film_20261102'
sys.path.insert(0,str(R/'work/plot_deps'))
import numpy as np
from PIL import Image,ImageDraw,ImageFont
from matplotlib.colors import SymLogNorm
from matplotlib import colormaps
PLAN=json.loads((H/'film_plan.json').read_text('utf8'))
ASSETS=json.loads((H/'layout_assets.json').read_text('utf8'))
GRAPH=json.loads((O/'graph_layout.json').read_text('utf8'))
with np.load(O/'graph_cache.npz',allow_pickle=False) as z:TIMES=z['time_s'];DRAGS=z['drag_N']
FONT='C:/Windows/Fonts/malgun.ttf';FONTS={};LOGOS={}
COLORS=['#67e8f9','#fda4af','#fb923c','#a3e635','#c4b5fd','#38bdf8','#facc15','#34d399','#f472b6','#f87171']
def font(n):
    if n not in FONTS:FONTS[n]=ImageFont.truetype(FONT,n)
    return FONTS[n]
def text(draw,xy,label,size,fill='#f2f5f9',anchor='la'):
    draw.text(xy,label,font=font(size),fill=fill,anchor=anchor)
def fitted(draw,label,size,width):
    while size>14 and draw.textlength(label,font=font(size))>width:size-=1
    return font(size)
def logo(cid,height):
    key=(cid,height)
    if key not in LOGOS:
        c=next(c for c in ASSETS['cases'] if c['id']==cid)
        with Image.open(R/c['logo']) as src:im=src.convert('RGBA')
        im=im.crop(im.getchannel('A').getbbox());im.thumbnail((int(height*2.5),height),Image.Resampling.LANCZOS)
        LOGOS[key]=im
    return LOGOS[key]
def label(im,cid,label,x,y,width,size,height):
    lg=logo(cid,height);im.paste(lg,(x,y),lg);d=ImageDraw.Draw(im);f=fitted(d,label,size,width-lg.width-15)
    d.text((x+lg.width+15,y+height/2),label,font=f,fill='#e9f0f6',anchor='lm')
def stage(i):
    if i<363:return 'pressure',min(120,i//3),None
    if i<453:return 'Q',120,None
    hero=0 if i<723 else 1;return 'hero',i-(453 if hero==0 else 723),PLAN['hero_cases'][hero]
def tile(cid,i,mode):
    if mode=='Q':return H/'tilesQ'/(cid+'.png')
    if i==362:return H/'previews'/(cid+'_rear.png')
    return H/'tiles'/cid/('f'+str(i).zfill(3)+'.png')
def legend(im,box,vorticity=False,compact=False):
    x,y,w,h=box;d=ImageDraw.Draw(im);fs=19 if compact else 24
    heading='Q*=50 등가면 · 와도 |ω| [s⁻¹]' if vorticity else '표시용 SurfaceDrag [Pa]'
    text(d,(x,y),heading,fs,fill='#b8c9d9')
    rgb=np.round(colormaps['magma' if vorticity else 'jet'](np.linspace(0,1,w))[:,:3]*255).astype('u1')
    ramp=Image.fromarray(np.repeat(rgb[None,:,:],12,axis=0));im.paste(ramp,(x,y+31));ramp.close()
    norm=SymLogNorm(linthresh=5,linscale=.6,vmin=-120,vmax=120,base=10,clip=True)
    ticks=[10,30,100,300] if vorticity else [-120,-5,0,5,120]
    for v in ticks:
        a=(math.log10(v)-1)/(math.log10(300)-1) if vorticity else float(norm(v))
        xx=x+round(a*(w-1));d.line((xx,y+44,xx,y+48),fill='#8a9daf')
        text(d,(xx,y+50),str(v),15 if compact else 19,fill='#9fb3c7',anchor='mt')
def graph(im,t,box,hero=None):
    x,y,w,h=box;d=ImageDraw.Draw(im);small=im.width<1500;fs=18 if small else 22
    left=x+(70 if small else 95);top=y+12;pw=w-(90 if small else 125);ph=h-60;ymax=GRAPH['ymax_N']
    for value in range(0,int(ymax)+1,400):
        yy=top+ph*(1-value/ymax);d.line((left,yy,left+pw,yy),fill='#1d2b38',width=1);text(d,(left-10,yy),str(value),fs,fill='#8ca4bb',anchor='rm')
    for value in [1.08,1.62,2.16,2.70,3.24]:
        xx=left+pw*(value-1.08)/2.16;d.line((xx,top,xx,top+ph),fill='#142330');text(d,(xx,top+ph+10),f'{value:.2f}',fs,fill='#9fb3c7',anchor='mt')
    for j,c in enumerate(PLAN['cases']):
        if hero and c['id'] not in PLAN['hero_cases']:continue
        index=int(np.searchsorted(TIMES,t,side='right'));index=max(1,min(index,len(TIMES)))
        xx=left+pw*(TIMES[:index]-1.08)/2.16;yy=top+ph*(1-DRAGS[j,:index]/ymax)
        points=list(zip(xx,yy));color=COLORS[j]
        if len(points)>1:d.line(points,fill=color,width=3 if small else 4)
    current=left+pw*(min(3.24,max(1.08,t))-1.08)/2.16;d.line((current,top,current,top+ph),fill='#d2dce6',width=1)
    text(d,(left,top-8),'DRAG [N]',fs,fill='#c3d1e0',anchor='ls');text(d,(left+pw,top-8),'TIME [s]',fs,fill='#b9c9da',anchor='rs')
def curve_key(im,y,hero=None):
    d=ImageDraw.Draw(im);small=im.width<1500;cases=[(j,c) for j,c in enumerate(ASSETS['cases']) if not hero or c['id'] in PLAN['hero_cases']]
    if small and not hero:
        for k,(j,c) in enumerate(cases):
            x=30+(k%5)*210;yy=y+(k//5)*25;d.line((x,yy+10,x+20,yy+10),fill=COLORS[j],width=3);text(d,(x+27,yy),c['short'],15,fill='#b8c9d9')
    else:
        cw=im.width/len(cases)
        for k,(j,c) in enumerate(cases):
            x=25+k*cw;d.line((x,y+10,x+30,y+10),fill=COLORS[j],width=4);text(d,(x+40,y),c['short'],18 if not small else 20,fill='#b8c9d9')
def wide(i):
    mode,index,cid=stage(i);im=Image.new('RGB',(2560,1440),'black');d=ImageDraw.Draw(im)
    t=(6+.1*index)*.18 if mode=='pressure' else 3.24
    if mode!='hero':
        text(d,(35,10),'10 CARS  |  동일 조건, 같은 시간',48)
        text(d,(35,73),f'100 km/h  ·  t = {t:.3f} s  ·  '+('압력 항력 · 전면 → 후면' if mode=='pressure' else '최종 속도장 · Q-criterion / Vorticity'),27,fill='#abc0d3')
        for j,c in enumerate(PLAN['cases']):
            x=(j%5)*512;y=110+(j//5)*460
            with Image.open(tile(c['id'],i,mode)) as src:im.paste(src.resize((512,400),Image.Resampling.LANCZOS),(x,y+31))
            label(im,c['id'],c['label'],x+18,y+4,476,29,35)
            value=float(np.interp(t,TIMES,DRAGS[j]));text(d,(x+22,y+432),f'현재 {value:,.1f} N',25,fill=COLORS[j])
        legend(im,(90,1046,1050,70))
        if mode=='Q':legend(im,(1410,1046,1050,70),True)
        else:text(d,(1430,1068),'동일 축척 · 원본 차체 유지',26,fill='#94aec3')
    else:
        j=next(j for j,c in enumerate(PLAN['cases']) if c['id']==cid);c=PLAN['cases'][j];rank=1 if cid==PLAN['hero_cases'][0] else 10
        label(im,cid,c['label'],100,18,1500,62,76)
        text(d,(2450,36),f'{rank}위',66,fill=COLORS[j],anchor='rt')
        text(d,(100,128),f"평균항력 {c['mean_native_corrected_drag_N']:,.1f} N",40)
        text(d,(2450,128),f"항력 역비율 {c['inverse_mean_ratio']:.2f}x",44,fill=COLORS[j],anchor='rt')
        text(d,(640,227),'PRESSURE DRAG',29,fill='#aac0d4',anchor='mt');text(d,(1920,227),'Q-CRITERION / VORTICITY',29,fill='#aac0d4',anchor='mt')
        for col,key in enumerate(['pressure','Q']):
            path=H/'heroes'/cid/(key+str(index).zfill(3)+'.png')
            with Image.open(path) as src:im.paste(src.resize((1280,640),Image.Resampling.LANCZOS),(col*1280,265))
        text(d,(1280,948),'t = 3.240 s · 최종 유동장 회전뷰',26,fill='#8fa9be',anchor='mt')
        legend(im,(90,1006,1050,70));legend(im,(1410,1006,1050,70),True)
    graph(im,t,(20,1130,2520,235),cid);curve_key(im,1390,cid)
    text(d,(1280,1423),'동일 t*=6..18의 평균항력 기준 · 모델/격자 조건의 비교 · 표시용 압력 보간',16,fill='#667e93',anchor='mt')
    return im
def tall(i):
    mode,index,cid=stage(i);im=Image.new('RGB',(1080,1920),'black');d=ImageDraw.Draw(im)
    t=(6+.1*index)*.18 if mode=='pressure' else 3.24
    if mode!='hero':
        text(d,(540,18),'10 CARS · 같은 조건의 비교',42,anchor='mt')
        text(d,(540,79),f'100 km/h  ·  t = {t:.3f} s',28,fill='#abc0d3',anchor='mt')
        for j,c in enumerate(PLAN['cases']):
            x=(j%2)*540;y=125+(j//2)*280
            with Image.open(tile(c['id'],i,mode)) as src:im.paste(src.resize((298,233),Image.Resampling.LANCZOS),(x+121,y+33))
            label(im,c['id'],c['label'],x+30,y+4,480,24,30)
            text(d,(x+270,y+256),f'현재 {np.interp(t,TIMES,DRAGS[j]):,.1f} N',22,fill=COLORS[j],anchor='mt')
        if mode=='pressure':legend(im,(100,1533,880,70),False,True)
        else:legend(im,(100,1533,880,70),True,True)
        graph(im,t,(5,1625,1070,195));curve_key(im,1840)
    else:
        j=next(j for j,c in enumerate(PLAN['cases']) if c['id']==cid);c=PLAN['cases'][j];rank=1 if cid==PLAN['hero_cases'][0] else 10
        text(d,(540,15),f'{rank}위 · 평균항력',45,fill=COLORS[j],anchor='mt')
        label(im,cid,c['label'],100,84,880,48,55)
        text(d,(540,175),f"{c['mean_native_corrected_drag_N']:,.1f} N  |  {c['inverse_mean_ratio']:.2f}x",59,fill=COLORS[j],anchor='mt')
        for key,y in [('pressure',292),('Q',918)]:
            text(d,(540,y-36),'PRESSURE DRAG' if key=='pressure' else 'Q-CRITERION / VORTICITY',24,fill='#abc0d3',anchor='mt')
            with Image.open(H/'heroes'/cid/(key+str(index).zfill(3)+'.png')) as src:im.paste(src.resize((1080,540),Image.Resampling.LANCZOS),(0,y))
        legend(im,(50,1472,460,70),False,True);legend(im,(570,1472,460,70),True,True)
        graph(im,t,(5,1615,1070,200),cid);curve_key(im,1850,cid)
    text(d,(540,1901),'동일 모델·격자 조건 · 평균항력 기준',17,fill='#667e93',anchor='mt')
    return im
def main():
    start=int(sys.argv[1]);count=int(sys.argv[2]);assert 0<=start<993 and 1<=count<=30 and start+count<=993
    for i in range(start,start+count):
        for name,fn in [('wide',wide),('tall',tall)]:
            path=H/('frames_final_'+name)/('f'+str(i).zfill(4)+'.png');assert not path.exists()
            im=fn(i);im.save(path,compress_level=2);im.close()
    print(json.dumps({'first_frame':start,'count':count,'layouts':['wide','tall'],'CFD_calls':0}))
if __name__=='__main__':main()
