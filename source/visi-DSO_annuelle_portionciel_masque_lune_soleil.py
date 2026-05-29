#! /usr/bin/env python

import ephem
import datetime
import math
import csv
import re
from collections import defaultdict
import numpy as np
from scipy.interpolate import interp1d

# ce code donne pour chaque jour de l'année les objets visibles au-dessus d'un masque (masque.txt) et en-dessous d'une  
# hauteur max (hmax) pendant une durée minimale (min_duration); les résultats sont écrits dans deux fichiers 
# (par date et par objet)

# Définition des objets du ciel profond (DSO)

DSO = []

with open("DSO.txt", "r", encoding="utf-8") as f:
    for line in f:
        # Extraire les données à l'aide d'une expression régulière
        match = re.search(r'\("([^"]*)",\s*([\d.-]+),\s*([\d.-]+),\s*([\d.-]+),\s*([\d.-]+)\)', line)

        if match:
            nom, val1, val2, val3, val4 = match.groups()
            DSO.append((nom, float(val1), float(val2), float(val3), float(val4)))
# Fonction pour charger le masque d'horizon depuis un fichier
def load_horizon_mask(filename):
    azimuths = []
    min_altitudes = []
    
    with open(filename, 'r') as file:
        for line in file:
            if line.strip() and not line.startswith('#'):
                parts = line.strip().split(';')
                if len(parts) >= 2:
                    az = float(parts[0].strip())
                    alt = float(parts[1].strip())
                    azimuths.append(az)
                    min_altitudes.append(alt)
    
    # S'assurer que le masque couvre 360 degrés
    if azimuths[0] != 0 or azimuths[-1] != 359:
        # Ajouter des points pour fermer le cercle si nécessaire
        if azimuths[0] != 0:
            azimuths.insert(0, 0)
            min_altitudes.insert(0, min_altitudes[-1] if azimuths[-1] == 359 else min_altitudes[0])
        if azimuths[-1] != 359:
            azimuths.append(359)
            min_altitudes.append(min_altitudes[0])
    
    # Créer une fonction d'interpolation
    # Conversion en radians
    azimuths_rad = [math.radians(az) for az in azimuths]
    min_altitudes_rad = [math.radians(alt) for alt in min_altitudes]
    
    # Utiliser une interpolation circulaire pour l'azimut
    azimuths_ext = azimuths_rad + [azimuths_rad[0] + 2*math.pi]
    min_altitudes_ext = min_altitudes_rad + [min_altitudes_rad[0]]
    
    horizon_func = interp1d(azimuths_ext, min_altitudes_ext, kind='linear', bounds_error=False, fill_value=(min_altitudes_rad[0], min_altitudes_rad[0]))
    
    return horizon_func

# Fonction pour vérifier si un objet est au-dessus du masque d'horizon
def is_above_horizon(az, alt, horizon_func):
    # S'assurer que l'azimut est dans [0, 2π]
    az_rad = float(az)
    while az_rad < 0:
        az_rad += 2*math.pi
    while az_rad >= 2*math.pi:
        az_rad -= 2*math.pi
    
    # Obtenir l'altitude minimale pour cet azimut
    min_alt = horizon_func(az_rad)
    
    # Vérifier si l'altitude est supérieure au minimum requis
    return (float(alt) >= min_alt and float(alt) < hmax)

# Fonction pour convertir les coordonnées équatoriales en coordonnées horizontales
def eq_to_horiz(ra, dec, observer):
    body = ephem.FixedBody()
    body._ra = ephem.hours(ra)
    body._dec = ephem.degrees(dec)
    body.compute(observer)
    return body.az, body.alt

# Fonction pour calculer la séparation angulaire entre deux objets
def angular_separation(ra1, dec1, ra2, dec2):
    # Convertir en radians si nécessaire
    ra1, dec1, ra2, dec2 = float(ra1), float(dec1), float(ra2), float(dec2)
    
    # Calculer la séparation angulaire en utilisant la formule de la trigonométrie sphérique
    cos_separation = math.sin(dec1) * math.sin(dec2) + math.cos(dec1) * math.cos(dec2) * math.cos(ra1 - ra2)
    
    # Éviter les erreurs numériques
    cos_separation = min(1.0, max(-1.0, cos_separation))
    
    # Convertir en degrés
    separation_rad = math.acos(cos_separation)
    separation_deg = math.degrees(separation_rad)
    
    return separation_deg

# Configuration de l'observateur
def configure_observer(lat, lon, elevation):
    observer = ephem.Observer()
    observer.lat = ephem.degrees(lat)
    observer.lon = ephem.degrees(lon)
    observer.elevation = elevation  # en mètres
    return observer

# Fonction principale
def calculate_visibility(year, min_duration, horizon_func, lat, lon, elevation):
    observer = configure_observer(lat, lon, elevation)
    
    # Dictionnaires pour stocker les résultats
    results_by_day = defaultdict(list)
    results_by_object = defaultdict(list)
    
    # Créer un objet Lune
    moon = ephem.Moon()
    
    # Créer un objet Soleil
    sun = ephem.Sun()
    
    # Parcourir chaque jour de l'année
    start_date = datetime.datetime(year, 1, 1)
    for day_offset in range(366 if year % 4 == 0 else 365):  # Gestion des années bissextiles
        current_date = start_date + datetime.timedelta(days=day_offset)
        date_str = current_date.strftime("%Y-%m-%d")
        
        # Période d'observation: de 17h00 à 7h00 le lendemain
        obs_start = current_date.replace(hour=17, minute=0, second=0)
        obs_end = (current_date + datetime.timedelta(days=1)).replace(hour=7, minute=0, second=0)
        
        # Pour chaque objet DSO
        for dso_name, ra_hours, dec_degrees, size, magnitude in DSO:
            # Convertir RA de heures en radians et Dec de degrés en radians
            ra = ra_hours * math.pi / 12
            dec = dec_degrees * math.pi / 180
            
            # Créer un objet fixe pour le DSO
            dso = ephem.FixedBody()
            dso._ra = ephem.hours(ra)
            dso._dec = ephem.degrees(dec)
            
            # Vérifier la visibilité par intervalles de temps
            time_step = datetime.timedelta(minutes=5)  # Vérification toutes les 5 minutes
            current_time = obs_start
            
            # Variables pour suivre les périodes de visibilité
            visibility_periods = []
            current_period_start = None
            
            # Variables pour suivre les distances à la Lune
            moon_distances = []
            
            # Variables pour suivre les temps où le soleil est sous 5 degrés
            valid_times = []
            valid_coords = []  # Pour stocker les coordonnées (az, alt) correspondant aux temps valides
            
            while current_time <= obs_end:
                observer.date = ephem.Date(current_time)
                
                # Calculer les coordonnées de l'objet
                dso.compute(observer)
                az, alt = dso.az, dso.alt
                
                # Calculer les coordonnées de la Lune
                moon.compute(observer)
                
                # Calculer les coordonnées du Soleil
                sun.compute(observer)
                
                # Calculer la distance angulaire entre l'objet et la Lune
                moon_distance = angular_separation(dso.ra, dso.dec, moon.ra, moon.dec)
                
                # Vérifier si le soleil est sous l'horizon (plus de 10 degrés)
                sun_altitude = float(sun.alt) * 180 / math.pi  # Conversion en degrés
                sun_below_limit = sun_altitude < -10
                
                # Vérifier si l'objet est au-dessus du masque d'horizon
                if is_above_horizon(az, alt, horizon_func):
                    if current_period_start is None:
                        current_period_start = current_time
                    
                    # Stocker la distance à la Lune
                    moon_distances.append(moon_distance)
                    
                    # Stocker si le soleil est assez bas à ce moment
                    if sun_below_limit:
                        valid_times.append(current_time)
                        valid_coords.append((az, alt))  # Stocker les coordonnées pour ce temps valide
                else:
                    if current_period_start is not None:
                        # Une période de visibilité vient de se terminer
                        # Ne l'enregistrer que s'il y a au moins un moment où le soleil est assez bas
                        if valid_times:
                            # Calculer la durée totale où le soleil est bas
                            valid_duration = sum((valid_times[i+1] - valid_times[i]).total_seconds() for i in range(len(valid_times)-1) if (valid_times[i+1] - valid_times[i]).total_seconds() <= 300) + 300
                            
                            if valid_duration / 60 >= min_duration and moon_distances:
                                min_moon_distance = min(moon_distances)
                                max_moon_distance = max(moon_distances)
                                
                                # Utiliser les coordonnées du premier et dernier temps valide
                                start_az, start_alt = valid_coords[0]
                                end_az, end_alt = valid_coords[-1]
                                
                                visibility_periods.append({
                                    'start_time': valid_times[0],
                                    'end_time': valid_times[-1],
                                    'duration': valid_duration / 60,
                                    'start_az': start_az,
                                    'start_alt': start_alt,
                                    'end_az': end_az,
                                    'end_alt': end_alt,
                                    'min_moon_distance': min_moon_distance,
                                    'max_moon_distance': max_moon_distance
                                })
                        
                        current_period_start = None
                        moon_distances = []
                        valid_times = []
                        valid_coords = []
                
                current_time += time_step
            
            # Vérifier la dernière période si elle est en cours à la fin de l'observation
            if current_period_start is not None and valid_times:
                # Calculer la durée totale où le soleil est bas
                valid_duration = sum((valid_times[i+1] - valid_times[i]).total_seconds() for i in range(len(valid_times)-1) if (valid_times[i+1] - valid_times[i]).total_seconds() <= 300) + 300
                
                if valid_duration / 60 >= min_duration and moon_distances:
                    # Utiliser les coordonnées du premier et dernier temps valide
                    start_az, start_alt = valid_coords[0]
                    end_az, end_alt = valid_coords[-1]
                    
                    min_moon_distance = min(moon_distances)
                    max_moon_distance = max(moon_distances)
                    
                    visibility_periods.append({
                        'start_time': valid_times[0],
                        'end_time': valid_times[-1],
                        'duration': valid_duration / 60,
                        'start_az': start_az,
                        'start_alt': start_alt,
                        'end_az': end_az,
                        'end_alt': end_alt,
                        'min_moon_distance': min_moon_distance,
                        'max_moon_distance': max_moon_distance
                    })
            
            # Ajouter les périodes de visibilité aux résultats
            for period in visibility_periods:
                result = {
                    'date': date_str,
                    'object': dso_name.strip(),
                    'size': size,
                    'magnitude': magnitude,
                    'start_time': period['start_time'].strftime("%H:%M"),
                    'end_time': period['end_time'].strftime("%H:%M"),
                    'duration': round(period['duration'], 1),
                    'start_az': math.degrees(float(period['start_az'])),
                    'start_alt': math.degrees(float(period['start_alt'])),
                    'end_az': math.degrees(float(period['end_az'])),
                    'end_alt': math.degrees(float(period['end_alt'])),
                    'min_moon_distance': round(period['min_moon_distance'], 1),
                    'max_moon_distance': round(period['max_moon_distance'], 1)
                }
                results_by_day[date_str].append(result)
                results_by_object[dso_name.strip()].append(result)
    
    # Écrire les résultats dans des fichiers
    write_results_by_day(results_by_day, f"resultats_par_jour_{year}.csv")
    write_results_by_object(results_by_object, f"resultats_par_objet_{year}.csv")

# Fonction pour écrire les résultats par jour
def write_results_by_day(results, filename):
    with open(filename, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Date', 'Objet', 'Taille (arcmin)', 'Magnitude', 'Heure début', 'Azimut début', 
                         'Hauteur début', 'Heure fin', 'Azimut fin', 'Hauteur fin', 'Durée (min)',
                         'Distance min Lune (°)', 'Distance max Lune (°)'])
        
        for date, objects in sorted(results.items()):
            for obj in sorted(objects, key=lambda x: x['start_time']):
                writer.writerow([
                    obj['date'],
                    obj['object'],
                    obj['size'],
                    obj['magnitude'],
                    obj['start_time'],
                    round(obj['start_az'], 1),
                    round(obj['start_alt'], 1),
                    obj['end_time'],
                    round(obj['end_az'], 1),
                    round(obj['end_alt'], 1),
                    obj['duration'],
                    obj['min_moon_distance'],
                    obj['max_moon_distance']
                ])

# Fonction pour écrire les résultats par objet
def write_results_by_object(results, filename):
    with open(filename, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Objet', 'Date', 'Taille (arcmin)', 'Magnitude', 'Heure début', 'Azimut début', 
                         'Hauteur début', 'Heure fin', 'Azimut fin', 'Hauteur fin', 'Durée (min)',
                         'Distance min Lune (°)', 'Distance max Lune (°)'])
        
        for obj_name, dates in sorted(results.items()):
            for obj in sorted(dates, key=lambda x: x['date'] + x['start_time']):
                writer.writerow([
                    obj['object'],
                    obj['date'],
                    obj['size'],
                    obj['magnitude'],
                    obj['start_time'],
                    round(obj['start_az'], 1),
                    round(obj['start_alt'], 1),
                    obj['end_time'],
                    round(obj['end_az'], 1),
                    round(obj['end_alt'], 1),
                    obj['duration'],
                    obj['min_moon_distance'],
                    obj['max_moon_distance']
                ])

# Fonction pour créer un fichier de masque d'exemple
def create_sample_mask_file(filename):
    with open(filename, 'w', newline='', encoding='utf-8') as f:
        f.write("# Format: azimut ; hauteur minimale\n")
        f.write("0 ; 90\n")
        f.write("90 ; 90\n")
        f.write("100 ; 30\n")
        f.write("260 ; 30\n")
        f.write("270 ; 90\n")
        f.write("359 ; 90\n")

# Exemple d'utilisation
if __name__ == "__main__":
    # Définir l'année
    year = 2026
    
    # Définir la durée minimale de visibilité en minutes
    min_duration = 90
    
    # Définir hauteur max
    hmax = math.radians(80)
 
    
    # Nom du fichier de masque d'horizon
    mask_filename = "masque.txt"
    
    # Créer un fichier de masque d'exemple si nécessaire
    # create_sample_mask_file(mask_filename)
    
    # Charger le masque d'horizon
    horizon_func = load_horizon_mask(mask_filename)
    
    # Définir la position de l'observateur (latitude, longitude, élévation)
    lat = math.radians(43.42)  # Garravet
    lon = math.radians(0.9)
    elevation = 200  # mètres
   
    # Calculer la visibilité
    calculate_visibility(year, min_duration, horizon_func, lat, lon, elevation)
    print(f"Calcul terminé. Résultats enregistrés dans les fichiers resultats_par_jour_{year}.csv et resultats_par_objet_{year}.csv")