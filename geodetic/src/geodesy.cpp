//----------------------------------------------------------------------------------------------------------------------
///
/// \file       geodesy.cpp
/// \brief      Земные эллипсоиды, геодезические задачи (персчеты координат)
/// \date       06.11.19 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#include <geodesy.h>

namespace SPML /// Библиотека СБПМ
{
namespace Geodesy /// Геодезические функции
{
//----------------------------------------------------------------------------------------------------------------------
CEllipsoid::CEllipsoid()
{
    a = 0.0;
    b = 0.0;
    f = 0.0;
    invf = 0.0;
}

CEllipsoid::CEllipsoid( double semiMajorAxis, double semiMinorAxis, double inverseFlattening, bool isInvfDef )
{
    a = semiMajorAxis;
    invf = inverseFlattening;
    if( isInvfDef && ( Compare::IsZero( inverseFlattening ) || std::isinf( inverseFlattening ) ) ) {
        b = semiMajorAxis;
        f = 0.0;
    } else if ( isInvfDef ) {
        b = ( 1.0 - ( 1.0 / inverseFlattening ) ) * semiMajorAxis;
        f = 1.0 / inverseFlattening;
    } else {
        b = semiMinorAxis;
        f = 1.0 / inverseFlattening;
    }
}

//----------------------------------------------------------------------------------------------------------------------
void GEOtoRAD( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    float latStart, float lonStart, float latEnd, float lonEnd, float &d, float &az, float &azEnd )
{
    // Параметры эллипсоида:
    double a = ellipsoid.A();
    double b = ellipsoid.B();
    double f = ellipsoid.F();

    // По умолчанию Радианы:
    double _latStart = static_cast<double>( latStart );
    double _lonStart = static_cast<double>( lonStart );
    double _latEnd = static_cast<double>( latEnd );
    double _lonEnd = static_cast<double>( lonEnd );

    // При необходимости переведем входные данные в Радианы:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _latStart *= Convert::DgToRdD;
            _lonStart *= Convert::DgToRdD;
            _latEnd *= Convert::DgToRdD;
            _lonEnd *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

    if( Compare::AreEqual( a, b ) ) { // При расчете на сфере используем упрощенные формулы
        // Azimuth
        double fact1, fact2, fact3;
        fact1 = cos( _latEnd ) * sin( _lonEnd - _lonStart );
        fact2 = cos( _latStart ) * sin( _latEnd );
        fact3 = sin( _latStart ) * cos( _latEnd ) * cos( _lonEnd - _lonStart );                
        az = static_cast<float>( Convert::AngleTo360( atan2( fact1, fact2 - fact3 ), Units::TAngleUnit::AU_Radian ) ); // [рад] - Прямой азимут в начальной точке

        // ReverseAzimuth
        fact1 = cos( _latStart ) * sin( _lonEnd - _lonStart );
        fact2 = cos( _latStart ) * sin( _latEnd ) * cos( _lonEnd - _lonStart );
        fact3 = sin( _latStart ) * cos( _latEnd );        
        azEnd = static_cast<float>( Convert::AngleTo360( ( atan2( fact1, fact2 - fact3 ) ), Units::TAngleUnit::AU_Radian ) ); // [рад] - Прямой азимут в конечной точке

        // Distance
        double temp1, temp2, temp3;
        temp1 = sin( _latStart ) * sin( _latEnd );
        temp2 = cos( _latStart ) * cos( _latEnd ) * cos( _lonEnd - _lonStart );
        temp3 = temp1 + temp2;        
        d = static_cast<float>( ( acos( temp3 ) * a ) ); // [м]
    } else { // Для эллипсоида используем формулы Винсента
        double L = _lonEnd - _lonStart;

        double U1 = atan( ( 1.0 - f ) * tan( _latStart ) );
        double U2 = atan( ( 1.0 - f ) * tan( _latEnd ) );

        double sinU1 = sin( U1 );
        double cosU1 = cos( U1 );
        double sinU2 = sin( U2 );
        double cosU2 = cos( U2 );

        // eq. 13
        double lambda = L;
        double lambda_new = 0;
        int iterLimit = 100;

        double sinSigma = 0;
        double cosSigma = 0;
        double sigma = 0;
        double sinAlpha = 0;
        double cosSqAlpha = 0;
        double cos2SigmaM = 0;
        double c = 0;
        double sinLambda = 0;
        double cosLambda = 0;

        do {
            sinLambda = sin( lambda );
            cosLambda = cos( lambda );

            // eq. 14
            sinSigma = sqrt( ( ( cosU2 * sinLambda ) * ( cosU2 * sinLambda ) +
                ( cosU1 * sinU2 - sinU1 * cosU2 * cosLambda ) * ( cosU1 * sinU2 - sinU1 * cosU2 * cosLambda ) ) );
            if( Compare::IsZero( sinSigma ) ) { // co-incident points
                d = 0;
                az = 0;
                azEnd = 0;
                return;
            }

            // eq. 15
            cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;

            // eq. 16
            sigma = atan2( sinSigma, cosSigma );

            // eq. 17    Careful!  sin2sigma might be almost 0!
            sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
            cosSqAlpha = 1 - sinAlpha * sinAlpha;

            // eq. 18    Careful!  cos2alpha might be almost 0!
            cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha;

            if( std::isnan( cos2SigmaM ) ) {
                cos2SigmaM = 0; // equatorial line: cosSqAlpha = 0
            }

            // eq. 10
            c = ( f / 16 ) * cosSqAlpha * ( 4 + f * ( 4 - 3 * cosSqAlpha ) );

            lambda_new = lambda;

            // eq. 11 (modified)
            lambda = L + ( 1 - c ) * f * sinAlpha *
                ( sigma + c * sinSigma * ( cos2SigmaM + c * cosSigma * ( -1 + 2 * cos2SigmaM * cos2SigmaM ) ) );

        } while( fabs( ( lambda - lambda_new ) / lambda ) > 1e-15 && --iterLimit > 0 ); // see how much improvement we got

        double uSq = cosSqAlpha * ( a * a - b * b ) / ( b * b );

        // eq. 3
        double A = 1 + uSq / 16384 * ( 4096 + uSq * ( -768 + uSq * ( 320 - 175 * uSq ) ) );

        // eq. 4
        double B = uSq / 1024 * ( 256 + uSq * ( -128 + uSq * ( 74 - 47 * uSq ) ) );

        // eq. 6
        double deltaSigma = B * sinSigma *
            ( cos2SigmaM + ( B / 4 ) * ( cosSigma * ( -1 + 2 * cos2SigmaM * cos2SigmaM ) -
            ( B / 6 ) * cos2SigmaM * ( -3 + 4 * sinSigma * sinSigma ) * ( -3 + 4 * cos2SigmaM * cos2SigmaM ) ) );

        // eq. 19        
        d = static_cast<float>( b * A * ( sigma - deltaSigma ) ); // [м]

        // eq. 20
        az = static_cast<float>( Convert::AngleTo360( atan2( ( cosU2 * sinLambda ),
            ( cosU1 * sinU2 - sinU1 * cosU2 * cosLambda ) ), Units::TAngleUnit::AU_Radian ) ); // Прямой азимут в начальной точке, [рад]

        // eq. 21
        azEnd = static_cast<float>( Convert::AngleTo360( atan2( ( cosU1 * sinLambda ),
            ( -sinU1 * cosU2 + cosU1 * sinU2 * cosLambda ) ), Units::TAngleUnit::AU_Radian ) ); // Прямой азимут в конечной точке, [рад]
    }
    // az, azEnd, d сейчас в радианах и метрах соответственно

    // Проверим, нужен ли перевод:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            az *= Convert::RdToDgF;
            azEnd *= Convert::RdToDgF;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            d *= 0.001f;
            break;
        }
        default:
            assert( false );
    }
    return;
}

void GEOtoRAD( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double latStart, double lonStart, double latEnd, double lonEnd, double &d, double &az, double &azEnd )
{
    // Параметры эллипсоида:
    double a = ellipsoid.A();
    double b = ellipsoid.B();
    double f = ellipsoid.F();

    // По умолчанию Радианы:
    double _latStart = latStart;
    double _lonStart = lonStart;
    double _latEnd = latEnd;
    double _lonEnd = lonEnd;

    // При необходимости переведем входные данные в Радианы:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _latStart *= Convert::DgToRdD;
            _lonStart *= Convert::DgToRdD;
            _latEnd *= Convert::DgToRdD;
            _lonEnd *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

    if( Compare::AreEqual( a, b ) ) { // При расчете на сфере используем упрощенные формулы
        // Azimuth
        double fact1, fact2, fact3;
        fact1 = cos( _latEnd ) * sin( _lonEnd - _lonStart );
        fact2 = cos( _latStart ) * sin( _latEnd );
        fact3 = sin( _latStart ) * cos( _latEnd ) * cos( _lonEnd - _lonStart );                
        az = Convert::AngleTo360( atan2( fact1, fact2 - fact3 ), Units::AU_Radian ); // [рад] - Прямой азимут в начальной точке

        // ReverseAzimuth
        fact1 = cos( _latStart ) * sin( _lonEnd - _lonStart );
        fact2 = cos( _latStart ) * sin( _latEnd ) * cos( _lonEnd - _lonStart );
        fact3 = sin( _latStart ) * cos( _latEnd );        
        azEnd = Convert::AngleTo360( ( atan2( fact1, fact2 - fact3 ) ), Units::AU_Radian ); // [рад] - Прямой азимут в конечной точке

        // Distance
        double temp1, temp2, temp3;
        temp1 = sin( _latStart ) * sin( _latEnd );
        temp2 = cos( _latStart ) * cos( _latEnd ) * cos( _lonEnd - _lonStart );
        temp3 = temp1 + temp2;        
        d = acos( temp3 ) * a ; // [м]        
    } else { // Для эллипсоида используем формулы Винсента
        double L = _lonEnd - _lonStart;

        double U1 = atan( ( 1.0 - f ) * tan( _latStart ) );
        double U2 = atan( ( 1.0 - f ) * tan( _latEnd ) );

        double sinU1 = sin( U1 );
        double cosU1 = cos( U1 );
        double sinU2 = sin( U2 );
        double cosU2 = cos( U2 );

        // eq. 13
        double lambda = L;
        double lambda_new = 0;
        int iterLimit = 100;

        double sinSigma = 0;
        double cosSigma = 0;
        double sigma = 0;
        double sinAlpha = 0;
        double cosSqAlpha = 0;
        double cos2SigmaM = 0;
        double c = 0;
        double sinLambda = 0;
        double cosLambda = 0;

        do {
            sinLambda = sin( lambda );
            cosLambda = cos( lambda );

            // eq. 14
            sinSigma = sqrt( ( ( cosU2 * sinLambda ) * ( cosU2 * sinLambda ) +
                ( cosU1 * sinU2 - sinU1 * cosU2 * cosLambda ) * ( cosU1 * sinU2 - sinU1 * cosU2 * cosLambda ) ) );
            if( Compare::IsZero( sinSigma ) ) { // co-incident points
                d = 0;
                az = 0;
                azEnd = 0;
                return;
            }

            // eq. 15
            cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;

            // eq. 16
            sigma = atan2( sinSigma, cosSigma );

            // eq. 17    Careful!  sin2sigma might be almost 0!
            sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
            cosSqAlpha = 1 - sinAlpha * sinAlpha;

            // eq. 18    Careful!  cos2alpha might be almost 0!
            cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha;

            if( std::isnan( cos2SigmaM ) ) {
                cos2SigmaM = 0; // equatorial line: cosSqAlpha = 0
            }

            // eq. 10
            c = ( f / 16 ) * cosSqAlpha * ( 4 + f * ( 4 - 3 * cosSqAlpha ) );

            lambda_new = lambda;

            // eq. 11 (modified)
            lambda = L + ( 1 - c ) * f * sinAlpha *
                ( sigma + c * sinSigma * ( cos2SigmaM + c * cosSigma * ( -1 + 2 * cos2SigmaM * cos2SigmaM ) ) );

        } while( fabs( ( lambda - lambda_new ) / lambda ) > 1e-15 && --iterLimit > 0 ); // see how much improvement we got

        double uSq = cosSqAlpha * ( a * a - b * b ) / ( b * b );

        // eq. 3
        double A = 1 + uSq / 16384 * ( 4096 + uSq * ( -768 + uSq * ( 320 - 175 * uSq ) ) );

        // eq. 4
        double B = uSq / 1024 * ( 256 + uSq * ( -128 + uSq * ( 74 - 47 * uSq ) ) );

        // eq. 6
        double deltaSigma = B * sinSigma *
            ( cos2SigmaM + ( B / 4 ) * ( cosSigma * ( -1 + 2 * cos2SigmaM * cos2SigmaM ) -
            ( B / 6 ) * cos2SigmaM * ( -3 + 4 * sinSigma * sinSigma ) * ( -3 + 4 * cos2SigmaM * cos2SigmaM ) ) );

        // eq. 19        
        d = b * A * ( sigma - deltaSigma ); // [m]

        // eq. 20
        az = Convert::AngleTo360( atan2( ( cosU2 * sinLambda ),
            ( cosU1 * sinU2 - sinU1 * cosU2 * cosLambda ) ), Units::TAngleUnit::AU_Radian ); // Прямой азимут в начальной точке, [рад]

        // eq. 21
        azEnd = Convert::AngleTo360( atan2( ( cosU1 * sinLambda ), ( -sinU1 * cosU2 + cosU1 * sinU2 * cosLambda ) ), Units::TAngleUnit::AU_Radian ); // Прямой азимут в конечной точке, [рад]
    }
    // az, azEnd, d сейчас в радианах и метрах соответственно

    // Проверим, нужен ли перевод:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            az *= Convert::RdToDgD;
            azEnd *= Convert::RdToDgD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            d *= 0.001;
            break;
        }
        default:
            assert( false );
    }
    return;
}
//----------------------------------------------------------------------------------------------------------------------
void RADtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    float latStart, float lonStart, float d, float az, float &latEnd, float &lonEnd, float &azEnd )
{        
    // Параметры эллипсоида:
    double a = ellipsoid.A();
    double b = ellipsoid.B();
    double f = ellipsoid.F();

    // по умолчанию Метры-Радианы:
    double _latStart = static_cast<double>( latStart ); // [рад]
    double _lonStart = static_cast<double>( lonStart ); // [рад]
    double _d = static_cast<double>( d );               // [м]
    double _az = static_cast<double>( az );             // [рад]

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _latStart *= Convert::DgToRdD;
            _lonStart *= Convert::DgToRdD;
            _az *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _d *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

    if( Compare::AreEqual( a, b ) ) { // При расчете на сфере используем упрощенные формулы
        _d = _d / a; // Нормирование
        // latitude
        double temp1, temp2, temp3;
        temp1 = sin( _latStart ) * cos( _d );
        temp2 = cos( _latStart ) * sin( _d ) * cos( _az );        
        latEnd = static_cast<float>( asin( temp1 + temp2 ) ); // [рад]

        // longitude
        temp1 = sin( _d ) * sin( _az );
        temp2 = cos( _latStart ) * cos( _d );
        temp3 = sin( _latStart ) * sin( _d ) * cos( _az );        
        lonEnd = static_cast<float>( _lonStart + atan2( temp1, temp2 - temp3 ) ); // [рад]

        // final bearing
        temp1 = cos( _latStart ) * sin( _az );
        temp2 = cos( _latStart ) * cos( _d ) * cos( _az );
        temp3 = sin( _latStart ) * sin( _d );        
        azEnd = static_cast<float>( Convert::AngleTo360( atan2( temp1, temp2 - temp3 ), Units::TAngleUnit::AU_Radian ) ); // [рад] - Прямой азимут в конечной точке
    } else { // Для эллипсоида используем формулы Винсента
        double cosAlpha1 = cos( _az );
        double sinAlpha1 = sin( _az );
        double s = _d; // distance [m]
        double tanU1 = ( 1.0 - f ) * tan( _latStart );
        double cosU1 = 1.0 / sqrt( ( 1.0 + tanU1 * tanU1 ) );
        double sinU1 = tanU1 * cosU1;

        // eq. 1
        double sigma1 = atan2( tanU1, cosAlpha1 );

        // eq. 2
        double sinAlpha = cosU1 * sinAlpha1;
        double cosSqAlpha = 1 - sinAlpha * sinAlpha;
        double uSq = cosSqAlpha * ( a * a - b * b ) / ( b * b );

        // eq. 3
        double A = 1 + ( uSq / 16384 ) * ( 4096 + uSq * ( -768 + uSq * ( 320 - 175 * uSq ) ) );

        // eq. 4
        double B = ( uSq / 1024 ) * ( 256 + uSq * ( -128 + uSq * ( 74 - 47 * uSq ) ) );

        // iterate until there is a negligible change in sigma
        double sOverbA = s / ( b * A );
        double sigma = sOverbA;
        double prevSigma = sOverbA;
        double cos2SigmaM = 0;
        double sinSigma = 0;
        double cosSigma = 0;
        double deltaSigma = 0;

        int iterations = 0;

        while( true ) {
            // eq. 5
            cos2SigmaM = cos( 2.0 * sigma1 + sigma );
            sinSigma = sin( sigma );
            cosSigma = cos( sigma );

            // eq. 6
            deltaSigma = B * sinSigma * ( cos2SigmaM +
                ( B / 4.0 ) * ( cosSigma * ( -1 + 2 * cos2SigmaM * cos2SigmaM ) -
                ( B / 6.0 ) * cos2SigmaM * ( -3 + 4 * sinSigma * sinSigma ) * ( -3 + 4 * cos2SigmaM * cos2SigmaM ) ) );

            // eq. 7
            sigma = sOverbA + deltaSigma;

            // break after converging to tolerance
            if( std::abs( sigma - prevSigma ) < 1e-15 || std::isnan( std::abs( sigma - prevSigma ) ) ) {
                break;
            }
            prevSigma = sigma;

            iterations++;
            if( iterations > 1000 ) {
                break;
            }
        }
        cos2SigmaM = cos( 2.0 * sigma1 + sigma );
        sinSigma = sin( sigma );
        cosSigma = cos( sigma );

        double tmp = sinU1 * sinSigma - cosU1 * cosSigma * cosAlpha1;

        // eq. 8
        latEnd = static_cast<float>( atan2( sinU1 * cosSigma + cosU1 * sinSigma * cosAlpha1,
            ( 1.0 - f ) * sqrt( ( sinAlpha * sinAlpha + tmp * tmp) ) ) ); // [рад]

        // eq. 9
        double lambda = atan2( ( sinSigma * sinAlpha1 ), ( cosU1 * cosSigma - sinU1 * sinSigma * cosAlpha1 ) );

        // eq. 10
        double c = ( f / 16 ) * cosSqAlpha * ( 4 + f * ( 4 - 3 * cosSqAlpha ) );

        // eq. 11
        double L = lambda - ( 1 - c ) * f * sinAlpha * ( sigma + c * sinSigma *
            ( cos2SigmaM + c * cosSigma * ( -1 + 2 * cos2SigmaM * cos2SigmaM ) ) );

        //double phi = ( _lonStart + L + 3 * PI ) % ( 2 * PI ) - PI;   // to -180.. 180 original!
//        lonEnd = static_cast<float>( ( _lonStart + L ) * RdToDgD );   // [град] my
        lonEnd = static_cast<float>( _lonStart + L );   // [рад] my

        // eq. 12
        double alpha2 = atan2( sinAlpha, -tmp ); // final direct bearing, if required        
        azEnd = static_cast<float>( Convert::AngleTo360( alpha2, Units::TAngleUnit::AU_Radian ) ); // Прямой азимут в конечной точке, [рад]
    }
    // latEnd, lonEnd, azEnd сейчас в радианах

    // Проверим, нужен ли перевод:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            latEnd *= Convert::RdToDgF;
            lonEnd *= Convert::RdToDgF;
            azEnd *= Convert::RdToDgF;
            break;
        }
        default:
            assert( false );
    }
    return;
}

void RADtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double latStart, double lonStart, double d, double az, double &latEnd, double &lonEnd, double &azEnd )
{
    // Параметры эллипсоида:
    double a = ellipsoid.A();
    double b = ellipsoid.B();
    double f = ellipsoid.F();

    // по умолчанию Метры-Радианы:
    double _latStart = latStart;    // [рад]
    double _lonStart = lonStart;    // [рад]
    double _d = d;                  // [м]
    double _az = az;                // [рад]

    // При необходимости переведем в Радианы-Метры:
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            _latStart *= Convert::DgToRdD;
            _lonStart *= Convert::DgToRdD;
            _az *= Convert::DgToRdD;
            break;
        }
        default:
            assert( false );
    }
    switch( rangeUnit ) {
        case( Units::TRangeUnit::RU_Meter ): break; // Уже переведено
        case( Units::TRangeUnit::RU_Kilometer):
        {
            _d *= 1000.0;
            break;
        }
        default:
            assert( false );
    }
    // Далее в математике используются углы в радианах и дальность в метрах, перевод в нужные единицы у конце

    if( Compare::AreEqual(a, b) ) { // При расчете на сфере используем упрощенные формулы
        _d = _d / a; // Нормирование
        // latitude
        double temp1, temp2, temp3;
        temp1 = sin( _latStart ) * cos( _d );
        temp2 = cos( _latStart ) * sin( _d ) * cos( _az );        
        latEnd = asin( temp1 + temp2 ); // [рад]

        // longitude
        temp1 = sin( _d ) * sin( _az );
        temp2 = cos( _latStart ) * cos( _d );
        temp3 = sin( _latStart ) * sin( _d ) * cos( _az );        
        lonEnd = _lonStart + atan2( temp1, temp2 - temp3 ); // [рад]

        // final bearing
        temp1 = cos( _latStart ) * sin( _az );
        temp2 = cos( _latStart ) * cos( _d ) * cos( _az );
        temp3 = sin( _latStart ) * sin( _d );        
        azEnd = Convert::AngleTo360( atan2( temp1, temp2 - temp3 ), Units::TAngleUnit::AU_Radian ); // [рад] - Прямой азимут в конечной точке
    } else { // Для эллипсоида используем формулы Винсента
        double cosAlpha1 = cos( _az );
        double sinAlpha1 = sin( _az );        
        double s = _d; // distance [m]
        double tanU1 = ( 1.0 - f ) * tan( _latStart );
        double cosU1 = 1.0 / sqrt( ( 1.0 + tanU1 * tanU1 ) );
        double sinU1 = tanU1 * cosU1;

        // eq. 1
        double sigma1 = atan2( tanU1, cosAlpha1 );

        // eq. 2
        double sinAlpha = cosU1 * sinAlpha1;
        double cosSqAlpha = 1 - sinAlpha * sinAlpha;
        double uSq = cosSqAlpha * ( a * a - b * b ) / ( b * b );

        // eq. 3
        double A = 1 + ( uSq / 16384 ) * ( 4096 + uSq * ( -768 + uSq * ( 320 - 175 * uSq ) ) );

        // eq. 4
        double B = ( uSq / 1024 ) * ( 256 + uSq * ( -128 + uSq * ( 74 - 47 * uSq ) ) );

        // iterate until there is a negligible change in sigma
        double sOverbA = s / ( b * A );
        double sigma = sOverbA;
        double prevSigma = sOverbA;
        double cos2SigmaM = 0;
        double sinSigma = 0;
        double cosSigma = 0;
        double deltaSigma = 0;

        int iterations = 0;

        while( true ) {
            // eq. 5
            cos2SigmaM = cos( 2.0 * sigma1 + sigma );
            sinSigma = sin( sigma );
            cosSigma = cos( sigma );

            // eq. 6
            deltaSigma = B * sinSigma * ( cos2SigmaM +
                ( B / 4.0 ) * ( cosSigma * ( -1 + 2 * cos2SigmaM * cos2SigmaM ) -
                ( B / 6.0 ) * cos2SigmaM * ( -3 + 4 * sinSigma * sinSigma ) * ( -3 + 4 * cos2SigmaM * cos2SigmaM ) ) );

            // eq. 7
            sigma = sOverbA + deltaSigma;

            // break after converging to tolerance
            if( std::abs( sigma - prevSigma ) < 1e-15 || std::isnan( std::abs( sigma - prevSigma ) ) ) {
                break;
            }
            prevSigma = sigma;

            iterations++;
            if( iterations > 1000 ) {
                break;
            }
        }
        cos2SigmaM = cos( 2.0 * sigma1 + sigma );
        sinSigma = sin( sigma );
        cosSigma = cos( sigma );

        double tmp = sinU1 * sinSigma - cosU1 * cosSigma * cosAlpha1;

        // eq. 8
        latEnd = atan2( sinU1 * cosSigma + cosU1 * sinSigma * cosAlpha1,
            ( 1.0 - f ) * sqrt( ( sinAlpha * sinAlpha + tmp * tmp) ) ); // [рад]

        // eq. 9
        double lambda = atan2( ( sinSigma * sinAlpha1 ), ( cosU1 * cosSigma - sinU1 * sinSigma * cosAlpha1 ) );

        // eq. 10
        double c = ( f / 16 ) * cosSqAlpha * ( 4 + f * ( 4 - 3 * cosSqAlpha ) );

        // eq. 11
        double L = lambda - ( 1 - c ) * f * sinAlpha * ( sigma + c * sinSigma *
            ( cos2SigmaM + c * cosSigma * ( -1 + 2 * cos2SigmaM * cos2SigmaM ) ) );

        //double phi = ( _lonStart + L + 3 * PI ) % ( 2 * PI ) - PI;   // to -180.. 180 original!
        //lonEnd = ( _lonStart + L ) * RdToDgD;   // [град] my
        lonEnd = _lonStart + L;   // [рад] my

        // eq. 12
        double alpha2 = atan2( sinAlpha, -tmp ); // final bearing, if required        
        azEnd = Convert::AngleTo360( alpha2, Units::TAngleUnit::AU_Radian ); // Прямой азимут в конечной точке, [рад]
    }
    // latEnd, lonEnd, azEnd сейчас в радианах

    // Проверим, нужен ли перевод:    
    switch( angleUnit ) {
        case( Units::TAngleUnit::AU_Radian ): break; // Уже переведено
        case( Units::TAngleUnit::AU_Degree ):
        {
            latEnd *= Convert::RdToDgD;
            lonEnd *= Convert::RdToDgD;
            azEnd *= Convert::RdToDgD;
            break;
        }
        default:
            assert( false );
    }
    return;
}

}
}
/// \}
