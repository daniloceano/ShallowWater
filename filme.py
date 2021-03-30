# -*- coding: utf-8 -*-

#IMPORTANDO PACOTES IMPORTANTES
from glob import glob
import imageio

def anim(path):
    filenames = glob(path+'*.png')
    ff=[]
    for i in range(len(filenames)):
    	ff.append(int(filenames[i].split('/')[-1].split('.')[0]))
    
    df = sorted(ff,key=int)
    
    filenames =[]
    for i in range(len(df)):
    	filenames.append(path+str(df[i])+'.png')
    
    with imageio.get_writer(path+'/leapfrogrotdiv.mp4',fps=5) as writer:
        for filename in filenames:
            image = imageio.imread(filename)
            writer.append_data(image)


