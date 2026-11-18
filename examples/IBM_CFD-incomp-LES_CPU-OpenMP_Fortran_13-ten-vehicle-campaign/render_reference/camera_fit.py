import itertools,json,math,sys
from pathlib import Path
H=Path(__file__).resolve().parent;R=H.parents[1]
sys.path.insert(0,str(R/'work/plot_deps'))
import numpy as np
def main():
    p=json.loads((H/'asset_plan.json').read_text('utf8'));boxes=[]
    for cid in ['bmw-e36-m3','land-rover-defender-110']:
        source=R/'outputs/car10_sns_smooth_20261102'/(cid+'_display_triangles.bin')
        report=json.loads((R/'outputs/car10_sns_smooth_20261102'/(cid+'_report.json')).read_text('utf8'))
        dtype=np.dtype([('parent','<i4'),('corner','<f8',(3,4))]);a=np.memmap(source,dtype=dtype,mode='r',offset=12,shape=(report['display_faces'],))
        lo=a['corner'][:,:,:3].min(axis=(0,1));hi=a['corner'][:,:,:3].max(axis=(0,1))
        boxes.append(np.asarray(list(itertools.product(*zip(lo,hi)))))
    scales={}
    for name,center,cap in [('pressure',np.array([0,1.15,0]),350),('Q',np.array([1.8,1.15,0]),220)]:
        allowed=float(cap)
        for azimuth in np.linspace(-140,400,541):
            el,az=np.radians([20,azimuth]);eye=np.array([np.cos(el)*np.cos(az),np.sin(el),np.cos(el)*np.sin(az)]);right=np.cross([0,1,0],eye);right/=np.linalg.norm(right);up=np.cross(eye,right)
            for corners in boxes:
                v=corners-center;x=v@right;y=v@up
                allowed=min(allowed,(1280-45)/max(abs(x)),(704-45)/max(y),(576-45)/max(-y))
        scales[name]=10*math.floor(allowed/10)
    with (R/'outputs/car10_final_film_20261102/camera_fit.json').open('x',encoding='utf8') as f:json.dump({'shared_hero_ppm_native':scales,'output_native_size':[2560,1280],'final_image_size':[1920,960],'same_scale_BMW_Defender':True,'all_bbox_corners_541_azimuths_fit_with_45px_margin':True,'actual_car_positions_modified':False},f,indent=2)
    print(json.dumps(scales))
if __name__=='__main__':main()
