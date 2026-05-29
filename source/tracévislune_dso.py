import matplotlib.pyplot as plt
import numpy as np
import re
import tkinter as tk
import pandas as pd
from astropy.time import Time
from astropy.coordinates import SkyCoord, AltAz, EarthLocation, solar_system_ephemeris, get_body
from astropy import units as u
from datetime import datetime, timedelta
from tkinter import ttk
from tkcalendar import Calendar
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from tkinter import filedialog
from matplotlib.path import Path
from matplotlib.patches import PathPatch

# ce programme trace les trajectoire de la lune et d'un DSO pour une date donnée et peut aussi prendre en compte un masque (masque.txt)

# Définition de quelques objets du ciel profond (nom, ascension droite, déclinaison)
DSO = []

with open("DSO_full.txt", "r", encoding="utf-8") as f:
    for line in f:
        # Extraire les données à l'aide d'une expression régulière
        match = re.search(r'\("([^"]*)",\s*([\d.-]+),\s*([\d.-]+),\s*([\d.-]+),\s*([\d.-]+)\)', line)

        if match:
            nom, val1, val2, val3, val4 = match.groups()
            DSO.append((nom, float(val1), float(val2), float(val3), float(val4)))
DEEP_SKY_OBJECTS = {
    dso[0].strip(): SkyCoord(ra=dso[1]*u.hour, dec=dso[2]*u.deg) 
    for dso in DSO
}

def interpolate_mask(mask_az, mask_alt):
    """
    Interpole les points du masque pour avoir une valeur par degré d'azimut
    
    Args:
        mask_az: tableau des azimuts du masque
        mask_alt: tableau des altitudes du masque
        
    Returns:
        az_interp: tableau des azimuts interpolés (un par degré)
        alt_interp: tableau des altitudes interpolées
    """
    # Trier les points par azimut croissant
    sort_idx = np.argsort(mask_az)
    mask_az = mask_az[sort_idx]
    mask_alt = mask_alt[sort_idx]
    
    # Créer un tableau d'azimuts avec un point par degré
    az_interp = np.arange(np.floor(min(mask_az)), np.ceil(max(mask_az)) + 1)
    
    # Interpoler les altitudes
    alt_interp = np.interp(az_interp, mask_az, mask_alt)
    
    return az_interp, alt_interp

def normalize_azimuth(az_array, center_az):
    """
    Normalise les azimuts autour d'un centre pour éviter la discontinuité 0-360
    """
    normalized = np.where(az_array - center_az > 180, az_array - 360, az_array)
    normalized = np.where(normalized - center_az < -180, normalized + 360, normalized)
    return normalized

def calculate_trajectories(date, start_hour, end_hour, observer_lat, observer_lon, deep_sky_object):
    """
    Calcule les trajectoires de la Lune et d'un objet du ciel profond
    """
    solar_system_ephemeris.set('builtin')
    location = EarthLocation(lat=observer_lat*u.deg, lon=observer_lon*u.deg)
    
    start_time = datetime.combine(date, datetime.strptime(start_hour, '%H:%M').time())
    end_time = datetime.combine(date, datetime.strptime(end_hour, '%H:%M').time())
    
    if end_time < start_time:
        end_time += timedelta(days=1)
    
    times = []
    current_time = start_time
    while current_time <= end_time:
        times.append(current_time)
        current_time += timedelta(minutes=5)
    
    times = Time(times)
    altaz_frame = AltAz(obstime=times, location=location)
    
    moon_coords = get_body('moon', times)
    moon_altaz = moon_coords.transform_to(altaz_frame)
    
    deep_sky_altaz = DEEP_SKY_OBJECTS[deep_sky_object].transform_to(altaz_frame)
    
    return moon_altaz, deep_sky_altaz, times

class TrajectoryApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Trajectoires célestes")
        
        # Paramètres par défaut
        self.observer_lat = 43.42  # Garravet
        self.observer_lon = 0.9
        
        # Création d'un conteneur principal avec poids
        self.root.grid_columnconfigure(0, weight=0)  # Contrôles (largeur fixe)
        self.root.grid_columnconfigure(1, weight=1)  # Graphique (extensible)
        
        # Frame pour les contrôles (à gauche)
        control_frame = ttk.Frame(root, padding="10")
        control_frame.grid(row=0, column=0, sticky="n")
        
        # Sélection de la date
        ttk.Label(control_frame, text="Date:").grid(row=0, column=0, sticky=tk.W)
        self.cal = Calendar(control_frame, selectmode='day', date_pattern='dd/mm/yyyy',
                          width=300, height=200)
        self.cal.grid(row=1, column=0, columnspan=2, pady=5)
        
        # Sélection des heures
        time_frame = ttk.Frame(control_frame)
        time_frame.grid(row=2, column=0, columnspan=2, pady=5)
        
        ttk.Label(time_frame, text="Heure début:").grid(row=0, column=0)
        self.start_hour = ttk.Entry(time_frame, width=8)
        self.start_hour.insert(0, "22:00")
        self.start_hour.grid(row=0, column=1, padx=5)
        
        ttk.Label(time_frame, text="Heure fin:").grid(row=0, column=2)
        self.end_hour = ttk.Entry(time_frame, width=8)
        self.end_hour.insert(0, "06:00")
        self.end_hour.grid(row=0, column=3, padx=5)
        
        # Sélection de l'objet
        ttk.Label(control_frame, text="Objet:").grid(row=3, column=0, sticky=tk.W, pady=5)
        self.object_var = tk.StringVar()
        self.object_combo = ttk.Combobox(control_frame, textvariable=self.object_var, 
                                       values=list(DEEP_SKY_OBJECTS.keys()), width=30)
        self.object_combo.grid(row=3, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=5)
        self.object_combo.set(list(DEEP_SKY_OBJECTS.keys())[0])
        
        # Boutons
        ttk.Button(control_frame, text="Calculer les trajectoires", 
                  command=self.plot_trajectories).grid(row=4, column=0, columnspan=2, pady=10)
        ttk.Button(control_frame, text="Sauvegarder PNG", 
                  command=self.save_plot).grid(row=5, column=0, columnspan=2, pady=5)
        ttk.Button(control_frame, text="Charger masque", 
                  command=self.load_mask).grid(row=6, column=0, columnspan=2, pady=5)
        
        self.mask_points = None
        
        # Frame pour le graphique (à droite)
        self.figure_frame = ttk.Frame(root)
        self.figure_frame.grid(row=0, column=1, sticky="nsew", padx=10, pady=10)
        
        # Configuration de la taille minimale de la fenêtre
        root.update_idletasks()
        min_width = control_frame.winfo_reqwidth() + 800  # Largeur minimale pour le graphique
        min_height = max(control_frame.winfo_reqheight(), 600)
        root.minsize(min_width, min_height)
        
    def load_mask(self):
        filename = filedialog.askopenfilename(filetypes=[("CSV files","*.csv"), ("Text files","*.txt")])
        if filename:
            mask_data = pd.read_csv(filename, sep=';', header=None)
            self.mask_points = mask_data.values
            self.plot_trajectories()

    def plot_trajectories(self):
        # Récupération des paramètres
        date_str = self.cal.get_date()
        date = datetime.strptime(date_str, '%d/%m/%Y').date()
        start_hour = self.start_hour.get()
        end_hour = self.end_hour.get()
        deep_sky_object = self.object_var.get()
        
        # Calcul des trajectoires
        moon_altaz, deep_sky_altaz, times = calculate_trajectories(
            date, start_hour, end_hour, self.observer_lat, self.observer_lon, deep_sky_object
        )
        
        # Création d'une nouvelle figure
        plt.close('all')
        fig = plt.figure(figsize=(8, 8))
        
        # Calcul du centre moyen des azimuts
        all_az = np.concatenate([moon_altaz.az.deg, deep_sky_altaz.az.deg])
        center_az = np.mean(all_az)
        
        # Normalisation des azimuts
        moon_az = normalize_azimuth(moon_altaz.az.deg, center_az)
        deep_sky_az = normalize_azimuth(deep_sky_altaz.az.deg, center_az)
        
                # Calcul des limites
        min_az = min(np.min(moon_az), np.min(deep_sky_az))
        max_az = max(np.max(moon_az), np.max(deep_sky_az))
        az_range = max_az - min_az
        az_margin = az_range * 0.01
        
        # Tracé
        plt.plot(moon_az, moon_altaz.alt.deg, 'y-', label='Lune')
        plt.plot(deep_sky_az, deep_sky_altaz.alt.deg, 'b-', label=deep_sky_object)
        

# Ajout du masque si présent
        if self.mask_points is not None:
            # D'abord interpoler les points du masque original
            az_interp, alt_interp = interpolate_mask(self.mask_points[:, 0], self.mask_points[:, 1])
            
            # Gérer la discontinuité 0°-360° pour le masque
            mask_az_normalized = []
            mask_alt_filtered = []
            
            for az, alt in zip(az_interp, alt_interp):
                # Normaliser chaque point d'azimut individuellement
                az_norm = normalize_azimuth(np.array([az]), center_az)[0]
                
                # Ne garder que les points dans la plage d'affichage
                if min_az - az_margin <= az_norm <= max_az + az_margin:
                    mask_az_normalized.append(az_norm)
                    mask_alt_filtered.append(alt)
            
            # Si nous avons des points de masque à afficher
            if mask_az_normalized:
                # Trier les points par azimut croissant
                mask_points_sorted = sorted(zip(mask_az_normalized, mask_alt_filtered))
                mask_az_sorted = [p[0] for p in mask_points_sorted]
                mask_alt_sorted = [p[1] for p in mask_points_sorted]
                
                # Créer le polygone du masque
                polygon_points = []
                
                # Points du masque de gauche à droite
                for az, alt in zip(mask_az_sorted, mask_alt_sorted):
                    polygon_points.append([az, alt])
                
                # Points du bas de droite à gauche pour fermer le polygone
                for az in reversed(mask_az_sorted):
                    polygon_points.append([az, -10])
                
                # Tracer le masque seulement s'il y a des points
                if len(polygon_points) > 2:
                    path = Path(polygon_points)
                    patch = PathPatch(path, facecolor='black', alpha=0.3)
                    plt.gca().add_patch(patch)        # Annotations temporelles
        total_points = len(times)
        step = int(total_points / (len(times) / 12))
        
        for i in range(0, total_points, step):
            time_str = times[i].datetime.strftime('%H:%M')
            plt.annotate(time_str, (moon_az[i], moon_altaz.alt.deg[i]),
                        xytext=(5, 5), textcoords='offset points',
                        fontsize=8, color='orange')
            plt.plot(moon_az[i], moon_altaz.alt.deg[i], 'yo', markersize=4)
            
            plt.annotate(time_str, (deep_sky_az[i], deep_sky_altaz.alt.deg[i]),
                        xytext=(-15, -5), textcoords='offset points',
                        fontsize=8, color='blue')
            plt.plot(deep_sky_az[i], deep_sky_altaz.alt.deg[i], 'bo', markersize=4)
        
        # Points cardinaux
        cardinal_points = [
            ('N', normalize_azimuth(np.array([0]), center_az)[0]),
            ('E', normalize_azimuth(np.array([90]), center_az)[0]),
            ('S', normalize_azimuth(np.array([180]), center_az)[0]),
            ('W', normalize_azimuth(np.array([270]), center_az)[0])
        ]
        
        for label, az in cardinal_points:
            if min_az - az_margin <= az <= max_az + az_margin:
                plt.text(az, -5, label, ha='center', fontsize=12, fontweight='bold')
        
        plt.xlim(min_az - az_margin, max_az + az_margin)
        plt.ylim(-10, 90)
        plt.xlabel('Azimut (degrés)')
        plt.ylabel('Altitude (degrés)')
        plt.title(f'Trajectoires de la Lune et de {deep_sky_object}\n'
                 f'Date: {date_str}\n'
                 f'De {start_hour} à {end_hour}')
        plt.grid(True)
        plt.legend()
        
        # Affichage dans la fenêtre tkinter
        for widget in self.figure_frame.winfo_children():
            widget.destroy()
        
        canvas = FigureCanvasTkAgg(fig, master=self.figure_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        
    def save_plot(self):
        date_str = self.cal.get_date().replace('/', '-')
        object_name = self.object_var.get().split()[0]
        filename = f"trajectoire_{object_name}_{date_str}.png"
        plt.savefig(filename, dpi=300, bbox_inches='tight')

if __name__ == "__main__":
    root = tk.Tk()
    app = TrajectoryApp(root)
    root.mainloop()