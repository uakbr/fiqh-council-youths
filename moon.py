from skyfield.api import load, Topos
from skyfield import almanac
from datetime import datetime, timedelta
import numpy as np
import math
import zipcodes  # For ZIP code to lat/long conversion

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

def get_location_from_zipcode(zip_code):
    """Convert ZIP code to latitude and longitude."""
    try:
        zip_data = zipcodes.matching(zip_code)
        if not zip_data:
            print(f"ZIP code {zip_code} not found. Using default location (Houston).")
            return 29.7604, -95.3698
        
        # Get the first match
        location = zip_data[0]
        lat = float(location['lat'])
        lon = float(location['long'])
        
        city = location['city']
        state = location['state']
        print(f"Location: {city}, {state}")
        
        return lat, lon
    except Exception as e:
        print(f"Error processing ZIP code: {e}")
        print("Using default location (Houston).")
        return 29.7604, -95.3698

def parse_coordinates(lat_str, lon_str):
    """Parse coordinates in various formats to decimal degrees."""
    # Handle latitude
    lat_str = lat_str.strip().upper()
    if 'N' in lat_str or 'S' in lat_str:
        lat_value = float(''.join([c for c in lat_str if c.isdigit() or c == '.' or c == '-']))
        if 'S' in lat_str:
            lat_value = -lat_value
    else:
        lat_value = float(lat_str)
        
    # Handle longitude
    lon_str = lon_str.strip().upper()
    if 'E' in lon_str or 'W' in lon_str:
        lon_value = float(''.join([c for c in lon_str if c.isdigit() or c == '.' or c == '-']))
        if 'W' in lon_str:
            lon_value = -lon_value
    else:
        lon_value = float(lon_str)
        
    return lat_value, lon_value

def get_user_location():
    """Prompt user for location using ZIP code or coordinates."""
    print("\nEnter your location:")
    print("1. US ZIP Code (e.g. 77429)")
    print("2. Coordinates (latitude and longitude)")
    
    choice = input("Choose option (1 or 2): ").strip()
    
    if choice == '1':
        zip_code = input("Enter ZIP code: ").strip()
        return get_location_from_zipcode(zip_code)
    elif choice == '2':
        print("Enter coordinates:")
        print("Examples: '29.7604 N', '29.7604', '29° 45' 37\" N'")
        lat_input = input("Latitude: ")
        lon_input = input("Longitude: ")
        
        try:
            return parse_coordinates(lat_input, lon_input)
        except ValueError:
            print("Invalid coordinates format. Using default location (Houston).")
            return 29.7604, -95.3698
    else:
        print(f"Invalid option: '{choice}'. Using default location (Houston).")
        return 29.7604, -95.3698

def moon_sighting_explanation(visible, details):
    """Provide a detailed explanation of moon sighting results tailored for Islamic calendar use."""
    explanation = []
    
    # Primary result explanation
    if visible:
        explanation.append("MOON SIGHTING RESULT: The lunar crescent SHOULD be visible.")
        explanation.append("This indicates the new Islamic month would likely begin.")
    else:
        explanation.append("MOON SIGHTING RESULT: The lunar crescent CANNOT be sighted.")
        explanation.append("This indicates the current Islamic month would likely continue.")
    
    # Detailed explanation based on parameters
    explanation.append("\nDetailed analysis:")
    
    # Altitude explanation
    if details['altitude'] < 0:
        explanation.append("- The moon is below the horizon at sunset (altitude: {:.1f}°)".format(details['altitude']))
        explanation.append("  This means the moon has already set before the sun, making sighting impossible.")
    else:
        if details['altitude'] < 5:
            explanation.append("- The moon is very low on the horizon (altitude: {:.1f}°)".format(details['altitude']))
            explanation.append("  Low altitude makes sighting difficult due to atmospheric extinction.")
        else:
            explanation.append("- The moon's altitude ({:.1f}°) is favorable for sighting.".format(details['altitude']))
    
    # Age explanation
    if details['moon_age_hours'] < 15:
        explanation.append("- The moon is very young ({:.1f} hours since new moon)".format(details['moon_age_hours']))
        explanation.append("  Crescents younger than 15-20 hours are extremely difficult to sight.")
    elif details['moon_age_hours'] < 24:
        explanation.append("- The moon is young ({:.1f} hours since new moon)".format(details['moon_age_hours']))
        explanation.append("  Young crescents require good atmospheric conditions to sight.")
    else:
        explanation.append("- The moon's age ({:.1f} hours) is favorable for sighting.".format(details['moon_age_hours']))
    
    # Elongation explanation
    if details['elongation'] < 7:
        explanation.append("- The angular separation between moon and sun ({:.1f}°) is below the Danjon limit (7°)".format(details['elongation']))
        explanation.append("  The crescent is physically impossible to see below this threshold.")
    elif details['elongation'] < 10:
        explanation.append("- The angular separation between moon and sun ({:.1f}°) is minimal".format(details['elongation']))
        explanation.append("  This small separation makes the crescent very thin and hard to sight.")
    else:
        explanation.append("- The angular separation ({:.1f}°) is favorable for sighting.".format(details['elongation']))
    
    # Illumination explanation
    illum_percent = details['illuminated_fraction'] * 100
    if illum_percent < 1:
        explanation.append("- The illuminated portion is extremely thin ({:.2f}%)".format(illum_percent))
        explanation.append("  Such a thin crescent requires exceptional viewing conditions.")
    elif illum_percent < 2:
        explanation.append("- The illuminated portion is very thin ({:.2f}%)".format(illum_percent))
        explanation.append("  A thin crescent requires clear skies and good eyesight.")
    else:
        explanation.append("- The illuminated portion ({:.2f}%) is favorable for sighting.".format(illum_percent))
    
    # Conclusion and recommendations
    explanation.append("\nConclusion:")
    if visible:
        if details['q_factor'] > 0.216:
            explanation.append("- The crescent should be easily visible to the naked eye")
            explanation.append("  Reliable sighting reports would be expected under normal conditions.")
        elif details['q_factor'] > -0.014:
            explanation.append("- The crescent may be visible under perfect atmospheric conditions")
            explanation.append("  Sighting reports should be carefully verified.")
        else:
            explanation.append("- The crescent may require optical aid to initially locate")
            explanation.append("  Naked-eye sighting reports should be treated with caution.")
    else:
        explanation.append("- The lunar crescent is not expected to be visible")
        explanation.append("  Any claimed sightings would be highly questionable.")
        # Add specific reasons
        if details['altitude'] < 0:
            explanation.append("  Primary reason: Moon below horizon at sunset")
        elif details['moon_age_hours'] < 15:
            explanation.append("  Primary reason: Moon too young ({:.1f} hours)".format(details['moon_age_hours']))
        elif details['elongation'] < 7:
            explanation.append("  Primary reason: Elongation below Danjon limit")
    
    # Join all explanations with newlines
    return "\n".join(explanation)

def main():
    """Main function to run the moon visibility check."""
    date_to_check = datetime.now()
    print(f"\nChecking lunar crescent visibility for {date_to_check.date()}")
    
    # Get user location
    lat, lon = get_user_location()
    
    # Display the location being used
    lat_dir = "N" if lat >= 0 else "S"
    lon_dir = "E" if lon >= 0 else "W"
    print(f"Using location: {abs(lat)}° {lat_dir}, {abs(lon)}° {lon_dir}")
    
    # Check moon visibility
    visible, sunset, details = check_moon_visibility(date_to_check, lat, lon)
    
    # Display results
    print(f"\nMoon visible on {date_to_check.date()}? {'Yes' if visible else 'No'}")
    if sunset:
        print(f"Sunset time (UTC): {sunset}")
        print(f"Details:")
        for key, value in details.items():
            print(f"  {key}: {value}")
    
    # Provide interpretations of the results
    print("\nInterpretation:")
    if details['altitude'] < 0:
        print("  - Moon is below the horizon at sunset")
    elif visible:
        if details['q_factor'] > 0.216:
            print("  - Crescent should be easily visible to the naked eye")
        elif details['q_factor'] > -0.014:
            print("  - Crescent may be visible in perfect atmospheric conditions")
        else:
            print("  - Crescent may require optical aid to initially locate")
    else:
        print("  - Lunar crescent not visible at this location/date")
        if details['moon_age_hours'] < 15:
            print("  - Moon is too young (less than 15 hours from new moon)")
        if details['elongation'] < 7:
            print("  - Elongation is below Danjon limit (7°)")

    # Ask if user wants to check another date or location
    while True:
        print("\nWould you like to check another date or location?")
        print("1. Try a different location")
        print("2. Try a different date")
        print("3. Exit")
        
        choice = input("Choose option (1, 2, or 3): ").strip()
        
        if choice == '1':
            lat, lon = get_user_location()
            visible, sunset, details = check_moon_visibility(date_to_check, lat, lon)
            main_display_results(date_to_check, lat, lon, visible, sunset, details)
        elif choice == '2':
            try:
                date_str = input("Enter date (YYYY-MM-DD): ").strip()
                new_date = datetime.strptime(date_str, "%Y-%m-%d")
                visible, sunset, details = check_moon_visibility(new_date, lat, lon)
                main_display_results(new_date, lat, lon, visible, sunset, details)
            except ValueError:
                print("Invalid date format. Please use YYYY-MM-DD.")
        else:
            break

def main_display_results(date, lat, lon, visible, sunset, details):
    """Display results for a specific date and location."""
    lat_dir = "N" if lat >= 0 else "S"
    lon_dir = "E" if lon >= 0 else "W"
    print(f"\nChecking lunar crescent visibility for {date.date()}")
    print(f"Location: {abs(lat)}° {lat_dir}, {abs(lon)}° {lon_dir}")
    
    print(f"\nMoon visible on {date.date()}? {'Yes' if visible else 'No'}")
    if sunset:
        print(f"Sunset time (UTC): {sunset}")
        print(f"Details:")
        for key, value in details.items():
            print(f"  {key}: {value}")
    
    # Provide comprehensive explanation
    print("\n" + "="*50)
    print(moon_sighting_explanation(visible, details))
    print("="*50)

if __name__ == "__main__":
    main()
