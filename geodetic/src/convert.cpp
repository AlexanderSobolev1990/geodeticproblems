//----------------------------------------------------------------------------------------------------------------------
///
/// \file       convert.cpp
/// \brief      Переводы единиц библиотеки СБПМ
/// \date       14.07.20 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#include <convert.h>

namespace SPML /// Библиотека СБПМ
{
namespace Convert /// Переводы единиц
{
//----------------------------------------------------------------------------------------------------------------------
float AngleTo360( float angle, const Units::TAngleUnit &au )
{    
    float _angle = angle;
    switch ( au ) {
        case Units::TAngleUnit::AU_Degree:
        {
            while( _angle >= 360.0f ) {
                _angle -= 360.0f;
            }
            while( _angle < 0.0f ) {
                _angle += 360.0f;
            }
            break;
        }
        case Units::TAngleUnit::AU_Radian:
        {
            while( _angle >= Consts::PI_2_F ) {
                _angle -= Consts::PI_2_F;
            }
            while( _angle < 0.0f ) {
                _angle += Consts::PI_2_F;
            }
            break;
        }
        default:
            assert( false );
    }
    return _angle;
}

double AngleTo360( double angle, const Units::TAngleUnit &au )
{    
    double _angle = angle;
    switch ( au ) {
        case Units::TAngleUnit::AU_Degree:
        {
            while( _angle >= 360.0 ) {
                _angle -= 360.0;
            }
            while( _angle < 0.0 ) {
                _angle += 360.0;
            }
            break;
        }
        case Units::TAngleUnit::AU_Radian:
        {
            while( _angle >= Consts::PI_2_F ) {
                _angle -= Consts::PI_2_F;
            }
            while( _angle < 0.0 ) {
                _angle += Consts::PI_2_F;
            }
            break;
        }
        default:
            assert( false );
    }
    return _angle;
}
//----------------------------------------------------------------------------------------------------------------------
float EpsToMP90( float angle, const Units::TAngleUnit &au )
{    
    float _angle = angle;
    switch ( au ) {
        case Units::TAngleUnit::AU_Degree:
        {
            while( _angle > 90.0f ) {
                _angle -= 180.0f;
            }
            while( _angle < -90.0f ) {
                _angle += 180.0f;
            }
            if( ( AngleTo360( std::abs( angle ), Units::TAngleUnit::AU_Degree ) > 90.0f ) &&
                ( AngleTo360( std::abs( angle ), Units::TAngleUnit::AU_Degree ) <= 270.0f ) )
            {
                 _angle *= -1.0f;
            }
            break;
        }
        case Units::TAngleUnit::AU_Radian:
        {
            while( _angle > Consts::PI_05_F ) {
                _angle -= Consts::PI_F;
            }
            while( _angle < -Consts::PI_05_F ) {
                _angle += Consts::PI_F;
            }
            if( ( AngleTo360( std::abs( angle ), Units::TAngleUnit::AU_Radian ) > Consts::PI_05_F ) &&
                ( AngleTo360( std::abs( angle ), Units::TAngleUnit::AU_Radian ) <= ( 3.0f * Consts::PI_05_F ) ) )
            {
                 _angle *= -1.0f;
            }
            break;
        }
        default:
            assert( false );
    }
    return _angle;
}

double EpsToMP90( double angle, const Units::TAngleUnit &au )
{    
    double _angle = angle;
    switch ( au ) {
        case Units::TAngleUnit::AU_Degree:
        {
            while( _angle > 90.0 ) {
                _angle -= 180.0;
            }
            while( _angle < -90.0 ) {
                _angle += 180.0;
            }
            if( ( AngleTo360( std::abs( angle ), Units::TAngleUnit::AU_Degree ) > 90.0 ) &&
                ( AngleTo360( std::abs( angle ), Units::TAngleUnit::AU_Degree ) <= 270.0 ) )
            {
                 _angle *= -1.0;
            }
            break;
        }
        case Units::TAngleUnit::AU_Radian:
        {
            while( _angle > Consts::PI_05_D ) {
                _angle -= Consts::PI_D;
            }
            while( _angle < -Consts::PI_05_D ) {
                _angle += Consts::PI_D;
            }
            if( ( AngleTo360( std::abs( angle ), Units::TAngleUnit::AU_Radian ) > Consts::PI_05_D ) &&
                ( AngleTo360( std::abs( angle ), Units::TAngleUnit::AU_Radian ) <= ( 3.0 * Consts::PI_05_D ) ) )
            {
                 _angle *= -1.0;
            }
            break;
        }
        default:
            assert( false );
    }
    return _angle;
}
//----------------------------------------------------------------------------------------------------------------------
void UnixTimeToHourMinSec(int rawtime, int &hour, int &min, int &sec, int &day, int &mon, int &year )
{
    std::time_t temp = rawtime;
    std::tm res;
    gmtime_r( &temp, &res );
    hour = ( res.tm_hour ) % 24;
    min = ( res.tm_min ) % 60;
    sec = ( res.tm_sec ) % 60;
    day = res.tm_mday;
    mon = ( res.tm_mon + 1 );
    year = ( res.tm_year + 1900 );
}
//----------------------------------------------------------------------------------------------------------------------
const std::string CurrentDateTimeToString() {
    time_t now = time( nullptr );
    struct tm tstruct;
    char buf[80];
    tstruct = *localtime(&now);
    // Visit http://en.cppreference.com/w/cpp/chrono/c/strftime
    // for more information about date/time format
    //strftime( buf, sizeof( buf ), "%Y-%m-%d.%X", &tstruct ); // original
    strftime( buf, sizeof( buf ), "%X %d-%m-%Y UTC%z", &tstruct ); // my
    return buf;
}

}
}
/// \}
