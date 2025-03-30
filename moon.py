from skyfield.api import load, Topos
from skyfield import almanac
from datetime import datetime, timedelta
import numpy as np
import math

# Load ephemeris data
eph = load('de421.bsp')
ts = load.timescale()

def check_moon_visibility(date, latitude, longitude):
    location = Topos(latitude, longitude)
    observer = eph['earth'] + location

    # Get sunset time
    t0 = ts.utc(date.year, date.month, date.day)
    t1 = ts.utc(date.year, date.month, date.day + 1)

    f = almanac.sunrise_sunset(eph, location)
    times, events = almanac.find_discrete(t0, t1, f)

    sunset_time = None
    for t, e in zip(times, events):
        if e == 0:  # 0 = sunset
            sunset_time = t
            break

    if sunset_time is None:
        # Handle polar regions where sun might not set
        return False, None

    # Check moon altitude at sunset
    moon = eph['moon']
    sun = eph['sun']
    
    # Moon position calculations
    moon_at_sunset = observer.at(sunset_time).observe(moon).apparent()
    moon_alt, moon_az, _ = moon_at_sunset.altaz()
    
    # Sun position at sunset
    sun_at_sunset = observer.at(sunset_time).observe(sun).apparent()
    sun_alt, sun_az, _ = sun_at_sunset.altaz()
    
    # Better elongation calculation using actual angular separation
    elongation = moon_at_sunset.separation_from(sun_at_sunset).degrees
    
    # Calculate moon age and illuminated fraction
    t_prev = ts.utc(date.year, date.month, date.day - 3)  # Look back 3 days to find new moon
    moon_phases = almanac.moon_phases(eph)
    phase_times, phase_events = almanac.find_discrete(t_prev, sunset_time, moon_phases)
    
    last_new_moon_time = None
    for t, e in zip(phase_times, phase_events):
        if e == 0:  # 0 = new moon
            last_new_moon_time = t
    
    moon_age_hours = 0
    if last_new_moon_time:
        moon_age_hours = (sunset_time.tt - last_new_moon_time.tt) * 24
    
    # Calculate illuminated fraction
    sun_to_moon = eph['moon'].at(sunset_time).position.km - eph['sun'].at(sunset_time).position.km
    earth_to_moon = eph['moon'].at(sunset_time).position.km - eph['earth'].at(sunset_time).position.km
    illum_fraction = (1 + np.dot(sun_to_moon, earth_to_moon) / 
                     (np.linalg.norm(sun_to_moon) * np.linalg.norm(earth_to_moon))) / 2
    
    # Calculate crescent width (approximation)
    crescent_width = illum_fraction * 2 * 1737.4  # Moon radius in km

    # Check moonset after sunset
    f_moon = almanac.risings_and_settings(eph, moon, location)
    moon_times, moon_events = almanac.find_discrete(t0, t1, f_moon)

    moonset_after_sunset = False
    for t, e in zip(moon_times, moon_events):
        if e == 1 and t > sunset_time:  # 1 = moonset
            moonset_after_sunset = True
            break

    # Calculate lag time (time between sunset and moonset)
    lag_time_minutes = 0
    for t, e in zip(moon_times, moon_events):
        if e == 1 and t > sunset_time:  # 1 = moonset
            lag_time_minutes = (t.tt - sunset_time.tt) * 24 * 60
            break
    
    # Yallop visibility criterion
    # q = (w - 0.00607 * a^2 + 0.1228 * a - 0.0702) / 10
    # where w is crescent width in arc minutes, a is geocentric elongation
    arc_min_per_km_at_moon_dist = 0.0145  # approximate conversion
    w = crescent_width * arc_min_per_km_at_moon_dist
    q = (w - 0.00607 * elongation**2 + 0.1228 * elongation - 0.0702) / 10
    
    # Danjon limit (minimum elongation based on moon age)
    min_elongation = max(7.0, 12.5 - 0.215 * moon_age_hours)
    
    # ARCL (ARCus Lucis) - arc of light along the lunar terminator
    arcl = elongation
    
    # Atmospheric extinction affects visibility
    extinction_factor = 1.0 / math.sin(math.radians(moon_alt.degrees + 2.0))
    
    # Determine visibility based on multiple criteria
    # Odeh criteria combined with traditional and Yallop's criteria
    visible = False
    
    # Basic requirements
    basic_requirements = (
        moonset_after_sunset and
        moon_alt.degrees >= 3 and  # Lower threshold, but we'll add more criteria
        elongation >= min_elongation
    )
    
    if basic_requirements:
        # Assign visibility based on criteria
        if q > 0.216:  # Easily visible to naked eye
            visible = True
        elif q > -0.014:  # Visible to naked eye in perfect conditions
            visible = True if extinction_factor < 2.5 and lag_time_minutes > 40 else False
        elif q > -0.160:  # May need optical aid to find, then naked eye
            visible = True if extinction_factor < 2.0 and moon_age_hours > 20 and lag_time_minutes > 50 else False
        elif q > -0.232:  # Optical aid needed to see
            visible = False  # Typically not visible to naked eye
        else:
            visible = False  # Not visible

    return visible, sunset_time.utc_datetime(), {
        'elongation': elongation,
        'moon_age_hours': moon_age_hours,
        'altitude': moon_alt.degrees,
        'illuminated_fraction': illum_fraction,
        'q_factor': q if 'q' in locals() else None,
        'lag_time_minutes': lag_time_minutes
    }

# Example: Houston on March 29, 2025
date_to_check = datetime(2025, 3, 29)
visible, sunset, details = check_moon_visibility(date_to_check, '29.7604 N', '95.3698 W')

print(f"Moon visible on {date_to_check.date()}? {visible}")
if sunset:
    print(f"Sunset time (UTC): {sunset}")
    print(f"Details:")
    for key, value in details.items():
        print(f"  {key}: {value}")
