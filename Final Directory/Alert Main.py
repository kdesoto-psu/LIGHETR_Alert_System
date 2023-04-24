import matplotlib.pyplot as plt
from astropy.time import Time
#import diagnose
#from diagnose import *
import json
import os
import pdb
from hop import Stream
from hop.io import StartPosition
import healpy as hp
from astropy.io import fits
import requests
import sys
sys.path.insert(0, '/Users/sky5265/Documents/GitHub/Astro_Code')
from Astro_useful_funcs import *
from Analysis_useful_funcs import *
from prob_obs import *
from get_galaxies import *
import get_galaxies
import get_LST
from twilio_caller_example import *
from twilio_texter_example import *
from testing_emailer import *
import time as TIme
import make_phaseii

recent_April = Time('2023-02-16T00:00:00.00')

def get_email_list(file_loc = 'contact_all_BNS.json'):
    f = open( file_loc , "rb" )
    jsonObject = json.load(f)
    return jsonObject['email_list']
    
def get_caller_list(file_loc = 'contact_all_BNS.json'):
    f = open( file_loc , "rb" )
    jsonObject = json.load(f)
    return jsonObject['caller_list']
    
def get_texter_list(file_loc = 'contact_all_BNS.json'):
    f = open( file_loc , "rb" )
    jsonObject = json.load(f)
    return jsonObject['texter_list']

def process_fits(recipients, fits_file, alert_message = None):
        
        alert_type = alert_message['alert_type']
        alert_time = alert_message['time_created']
        
        alert_time = Time(alert_time)
        
        if alert_time.jd < recent_April.jd:
            return
        
        #so we've found a time that we want to look at. I'll make a directory for this time.
        obs_time_dir = str(alert_time.mjd)+"/"
        mkdir(obs_time_dir)
        #I want to store the multi and single order fits here, along with the alert message as well.
        #write_to_file(obs_time_dir+"alert_message.txt", alert_message)
        


        #the fits file should be the multi
        
        
        #os.system('ligo-skymap-plot-observability '+str(fits_file)+' --site HET -o observability.png')
        #os.system('ligo-skymap-plot '+str(fits_file)+' -o '+str(alert_message['superevent_id'])+'.png')
        
        
        
        #skymap,header = hp.read_map(str(params['id'])+'.singleorder.fits', h=True, verbose=False)
        superevent_id = alert_message['superevent_id']
        fits_url = 'https://gracedb.ligo.org/api/superevents/'+str(superevent_id)+'/files/bayestar.multiorder.fits'
        
        
        #download the multi-order fits file from the fits_url filelocation
        url = fits_url
        multiorder_file_name = obs_time_dir+'multiorder_fits_'+str(superevent_id)+'.fits'
        req = requests.get(url)
        file = open(multiorder_file_name, 'wb')
        for chunk in req.iter_content(100000):
            file.write(chunk)
        file.close()
        
        #flatten the multi-order fits file into single-order
        singleorder_file_name = obs_time_dir+'flattened_multiorder_fits_'+superevent_id+'.fits'
        os.system('ligo-skymap-flatten '+str(multiorder_file_name)+' '+singleorder_file_name)
        
        #open the flattened fits file
        fits_file = fits.open(singleorder_file_name)
        
        #get the skymap and header of the flattened, single-order fits file
        skymap,header = hp.read_map(singleorder_file_name, h=True, verbose=False)
        
        
        #plot the sky-localization from the flattened, single-order fits file
        m1 = hp.read_map(singleorder_file_name)
        hp.mollview(m1, rot = (180, -10, 0))
        plt.savefig(obs_time_dir+'fits_plotted.png')
        plt.close()
        
        
        
        
        #skymap,header = hp.read_map(str(fits_url), h=True, verbose=False)
        #params['skymap_fits'] = fits_file
        

        # Print and save some values from the FITS header.
        header = dict(header)
        
        
        obs_time = Time(header['DATE-OBS'],format='isot',scale='utc')
        #time = Time.now()
        time = Time('2023-04-17T07:00:00.00')
        dist = str(header['DISTMEAN']) + ' +/- ' + str(header['DISTSTD'])
        header['id'] = superevent_id

        # Making a pie chart of the type of event for the email
        
        
        
        
        
        '''
        print("Alert message event keys: "+str(alert_message['event'].keys()))
        print("Alert message keys: "+str(alert_message.keys()))
        print("Alert message event properties keys: "+str(alert_message['event']['properties'].keys()))
        print("alert emssage event properties hasNs:" +str(alert_message['event']['properties']['HasNS']))
        print("alert emssage event properties hasRemnant:" +str(alert_message['event']['properties']['HasRemnant']))
        print("alert emssage event properties hasNs:" +str(alert_message['event']['classification']))
        
        
        The format of these alerts is given in this website:
        https://emfollow.docs.ligo.org/userguide/content.html
        
        So, the probabilities of BBH, BNS, and so on is located in:
        alert_message['event']['classification']['BNS']
        alert_message['event']['classification']['BBH']
        alert_message['event']['classification']['NSBH']
        alert_message['event']['classification']['Noise']
        
        it's unclear whether the key for noise will be either 'Noise' or 'Terrestrial', so I'll make the code flexible enough to use both
        
        '''
        noise_or_terrestrial = 'Terrestrial'
        if noise_or_terrestrial not in  alert_message['event']['classification'].keys():
            noise_or_terrestrial = 'Noise'
        event_poss = ['BBH', 'BNS', 'NSBH', noise_or_terrestrial]
        labels = []
        sizes = []
        
        for label in event_poss:
            val = float(alert_message['event']['classification'][label])
            if val >= 0:
                sizes.append(val)
                labels.append(label)
        
        
        #labels = ['%s (%.1f %%)'%(lab,pct) for lab,pct in zip(labels,sizes)]
        fig1, ax1 = plt.subplots()
        patches,texts = ax1.pie(sizes, startangle=90)
        ax1.legend(patches, labels, loc="best")
        ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
        plt.savefig(obs_time_dir+'piechart.png')

        

        prob, probfull, timetill90, m = prob_observable(skymap, header, time, savedir = obs_time_dir, plot=True)

        alert_message['skymap_fits'] = singleorder_file_name
        alert_message['skymap_array'] = m
        
        
        #if False:
        if timetill90 ==-99:
            print("HET can't observe the source.")
            
            return
        else:
            print("Source has a {:.1f}% chance of being observable now.".format(int(round(100 * prob))))
            print("Integrated probability over 24 hours (ignoring the sun) is {:.1f}%".format(int(round(100 * probfull))))
            print('{:.1f} hours till you can observe the 90 % prob region.'.format(timetill90))
            write_to_file(obs_time_dir+" observability.txt", "Source has a {:.1f}% chance of being observable now.".format(int(round(100 * prob))), append = False)
            write_to_file(obs_time_dir+" observability.txt", "Integrated probability over 24 hours (ignoring the sun) is {:.1f}%".format(int(round(100 * probfull))), append = True)
            write_to_file(obs_time_dir+" observability.txt", '{:.1f} hours till you can observe the 90 % prob region.'.format(timetill90), append = True)
            
            cattop, logptop, num_galaxies_visible_HET = get_galaxies.write_catalog(alert_message, 'GLADE', savedir = obs_time_dir)
            print("cattop: "+str(cattop))
            mincontour = get_LST.get_LST(savedir = obs_time_dir,targf = obs_time_dir+'galaxies%s_%s.dat'%('GLADE',alert_message['superevent_id']))
            
            
            
            contact_list_file_loc = 'contact_only_HET_BNS.json'
            
            #email(contact_list_file_loc = contact_list_file_loc, subject='LIGHETR Alert: NS Merger Detected', body = 'A Neutron Star Merger has been detected by LIGO.\n{:.1f} hours till you can observe the 90 % prob region.'.format(timetill90)+"\nSource has a {:.1f}% chance of being observable now.".format(int(round(100 * prob))), files_to_attach = [obs_time_dir+"HET_visibility_figure.png", obs_time_dir+"piechart.png"], people_to_contact = ['Karthik'])
            
            
            
            email_body = 'A Neutron Star Merger has been detected by LIGO.\n{:.1f} hours till you can observe the 90 % prob region.'.format(timetill90)+"\nSource has a {:.1f}% chance of being observable now.".format(int(round(100 * prob)))
            #email_body = 'Hi'
            
            email(contact_list_file_loc = contact_list_file_loc, subject='LIGHETR Alert: NS Merger Detected', body = email_body, files_to_attach = [], people_to_contact = ['Karthik', 'Karthik_PSU'])
            
            
            calling_dict = get_caller_list(contact_list_file_loc)
            texting_dict = get_texter_list(contact_list_file_loc)
            
            make_phaseii.make_phaseii(obs_time_dir+'LSTs_{}.out'.format(alert_message['superevent_id']))
            
            call_people(calling_dict = calling_dict, people_to_contact = ['Karthik'], message_to_say = 'NS Event Detected. Check email for information.')
            send_text_messages(reciever_dict = texting_dict, people_to_contact = ['Karthik'], message_to_send = 'NS Event Detected. Check email for information.')
            
            #send_notifications(params,timetill90)
            
            
            
###########Things start here####################



#stream_start_pos = 1600
stream_start_pos = StartPosition.EARLIEST
#print("Starting stream at "+str(stream_start_pos))
stream = Stream(start_at=stream_start_pos)

#stream = Stream()

num_messages = 0

print("Listening for Alerts from kafka")

with stream.open("kafka://kafka.scimma.org/igwn.gwalert", "r") as s:
    for message in s:
        print("Got an alert")
        event = message.content[0]['event']
        message_content = message.content[0]
        alert_type = message_content['alert_type']
        alert_time = message_content['time_created']
        
        alert_time = Time(alert_time)
        
        
        
        
        if alert_time.jd < recent_April.jd:
            num_messages+=1
            continue
        
        
        if event is not None:
            print("Event Keys: "+str(event.keys()))
            
        
        skymap = None
        if 'skymap' in message_content.keys():
            skymap = message_content['skymap']
        if event is not None and 'skymap' in event.keys():
            skymap = event['skymap']

        '''send that fits file into process_fits'''
        if skymap is not None and event is not None:
            print('Calling process_fits')
            process_fits(recipients = 'recipients.py', fits_file = skymap, alert_message = message_content)




#decide if it's visible to het

#if yes, send emails/text/phone calls to people

#create TSL for the files
