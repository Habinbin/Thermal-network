import math
from datetime import datetime
import constant as c


# 각도 단위 변경
d2r = math.pi/180
r2d = 180/math.pi

# 서울 위도, 경도
local_latitude = 37.478
local_longitude = 127.148
standard_longitude = 135

# 날짜 / 시간
year = 2024
month = 6
day = 22
local_hour = 12
local_min = 0

def equation_of_time(day_of_year): # Equation of Time: 진태양시와 평균태양시의 차이
    B = (day_of_year - 1) * 360/365 # B: 각도 단위

    EOT = 229.2 * (0.000075
               + 0.001868 * math.cos(d2r * B)
               - 0.032077 * math.sin(d2r * B)
               - 0.014615 * math.cos(d2r * 2 * B)
               - 0.04089 * math.sin(d2r * 2 * B)) 
    return EOT


def solar_position(year, month, day, local_hour, local_min, local_latitude, local_longitude, standard_longitude):
    # Equation of Time
    day_of_year = datetime(year, month, day).timetuple().tm_yday # 1월 1일부터 몇 번째 날인지
    EOT = equation_of_time(day_of_year)
    
    # Solar Time
    local_hour_decimal = (local_hour + local_min)*c.m2h # [h]
    delta_longitude = local_longitude - standard_longitude # [deg]
    solar_time_decimal = local_hour_decimal + (4 * delta_longitude + EOT) * c.m2h # [h]
    solar_time_min = (solar_time_decimal * 60) % 60 # [min]
    
    hour_angle = solar_time_decimal * 15 - 180 # [deg]
    
    # Solar Declination
    solar_declination = 23.45 * math.sin(d2r * 360 / 365 * (284 + day_of_year))
    
    # Solar Altitude, Solar Azimuth
    term_1 = math.cos(d2r * local_latitude) * math.cos(d2r * solar_declination) * math.cos(d2r * hour_angle) \
             + math.sin(d2r * local_latitude) * math.sin(d2r * solar_declination)
    solar_altitude = r2d * math.asin(term_1)
    
    
    term_2 = (math.sin(d2r * solar_altitude) * math.sin(d2r * local_latitude) - math.sin(d2r * solar_declination))\
             / (math.cos(d2r * solar_altitude) * math.cos(d2r * local_latitude))
    solar_azimuth = r2d * math.acos(term_2)
    
    # solar_zenith_angle = r2d * math.acos(math.sin(d2r * local_latitude) * math.sin(d2r * solar_declination) + math.cos(d2r * local_latitude) * math.cos(d2r * solar_declination) * math.cos(d2r * hour_angle))
    term_3 = (math.sin(d2r * solar_declination) * math.cos(d2r * local_latitude) - math.cos(d2r * hour_angle) * math.cos(d2r * solar_declination) * math.sin(d2r * local_latitude)) / math.sin(d2r * (90-solar_altitude))
    solar_azimuth_2 = r2d * math.acos(term_3)
    return solar_altitude, solar_azimuth