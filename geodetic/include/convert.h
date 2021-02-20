//----------------------------------------------------------------------------------------------------------------------
///
/// \file       convert.h
/// \brief      Переводы единиц библиотеки СБПМ
/// \date       27.07.20 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#ifndef SPML_CONVERT_H
#define SPML_CONVERT_H

#include <cmath>
#include <ctime>
#include <string>
#include <cassert>
#include <type_traits>

// SPML includes:
#include <consts.h>
#include <units.h>

namespace SPML /// Библиотека СБПМ
{
namespace Convert /// Переводы единиц
{
//----------------------------------------------------------------------------------------------------------------------
// Константы перевода радианов в градусы и наоборот
const float DgToRdF = static_cast<float>( asin( 1.0 ) / 90.0 ); ///< Перевод градусов в радианы (float) (путем умножения на данную константу)
const float RdToDgF = static_cast<float>( 90.0 / asin( 1.0 ) ); ///< Перевод радианов в градусы (float) (путем умножения на данную константу)
const double DgToRdD = asin( 1.0 ) / 90.0;                      ///< Перевод градусов в радианы (double) (путем умножения на данную константу)
const double RdToDgD = 90.0 / asin( 1.0 );                      ///< Перевод радианов в градусы (double) (путем умножения на данную константу)

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Приведение угла в [0,360) градусов или [0,2PI) радиан
/// \param[in] angle - приводимый угол в [град] или [рад] в зависимости от параметра.
/// \param[in] au    - выбор угловых единиц [град] или [рад]
/// \return Значение angle, приведенное в [0,360) градусов или [0,2PI) радиан
///
float AngleTo360( float angle, const Units::TAngleUnit &au );

///
/// \brief Приведение угла в [0,360) градусов или [0,2PI) радиан
/// \param[in] angle - приводимый угол в [град] или [рад] в зависимости от параметра.
/// \param[in] au    - выбор угловых единиц [град] или [рад]
/// \return Значение angle, приведенное в [0,360) градусов или [0,2PI) радиан
///
double AngleTo360( double angle, const Units::TAngleUnit &au );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Приведение угла места в [-90,90] градусов или [-PI/2, PI/2] радиан
/// \param[in] angle - приводимый угол в [град] или [рад] в зависимости от параметра.
/// \param[in] au    - выбор угловых единиц [град] или [рад]
/// \return Значение angle, приведенное в [-90,90] градусов или [-PI/2, PI/2] радиан
///
float EpsToMP90( float angle, const Units::TAngleUnit &au );

///
/// \brief Приведение угла места в [-90,90] градусов или [-PI/2, PI/2] радиан
/// \param[in] angle - приводимый угол в [град] или [рад] в зависимости от параметра.
/// \param[in] au    - выбор угловых единиц [град] или [рад]
/// \return Значение angle, приведенное в [-90,90] градусов или [-PI/2, PI/2] радиан
///
double EpsToMP90( double angle, const Units::TAngleUnit &au );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Перевод абсолютного азимута относительно севера в азимут относительно указанного направления
/// \param[in] absAz  - абсолютный азимут относительно севера
/// \param[in] origin - абсолютный азимут, относительно которого измеряется относительный
/// \param[in] au     - единицы измерения углов входных/выходных параметров
/// \return Азимут относительно указанного направления origin (положительный азимут 0..180 по часовой стрелке, отрицательный 0..-180 против часовой стрелки)
///
template <class T> T AbsAzToRelAz( T absAz, T origin, const Units::TAngleUnit &au )
{
    assert(
        ( std::is_same<T, float>::value ) ||
        ( std::is_same<T, double>::value )
    );
    T result = AngleTo360( absAz, au ) - origin ;
    return result;
}

///
/// \brief Перевод относительного азимута в абсолютный азимут относительно севера
/// \details Единицы измерения углов входных/выходных параметров согласно TAngleUnit
/// \param[in] relAz  - относительный азимут
/// \param[in] origin - абсолютный азимут, относительно которого измеряется relAz
/// \param[in] au     - единицы измерения углов входных/выходных параметров
/// \return Абсолютный азимут относительно севера
///
template <class T> T RelAzToAbsAz( T relAz, T origin, const Units::TAngleUnit &au )
{
    assert(
        ( std::is_same<T, float>::value ) ||
        ( std::is_same<T, double>::value )
    );
    T result = AngleTo360( ( relAz + origin ), au );
    return result;
}

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Перевод [дБ] в разы по мощности
/// \param[in] dB - децибелы
/// \return децибелы, перевденные в разы по мощности
///
template <class T> inline T dBtoTimesByP( T dB )
{
    assert(
        ( std::is_same<T, float>::value ) ||
        ( std::is_same<T, double>::value )
    );
    return ( std::pow( 10.0, dB / 10.0 ) );
}

///
/// \brief Перевод [дБ] в разы по напряжению
/// \param[in] dB - децибелы
/// \return децибелы, перевденные в разы по напряжению
///
template <class T> inline T dBtoTimesByU( T dB )
{
    assert(
        ( std::is_same<T, float>::value ) ||
        ( std::is_same<T, double>::value )
    );
    return ( std::pow( 10.0, dB / 20.0 ) );
}

static int dummy_int;
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Перевод целого числа секунд с 00:00:00 01.01.1970 в часы/минуты/секунды/день/месяц/год
/// \param[in] rawtime - число секунд с 00:00:00 01.01.1970
/// \param[out] hour - часы
/// \param[out] min - минуты
/// \param[out] sec - секунды
/// \param[out] day - день
/// \param[out] mon - месяц
/// \param[out] year - год
///
void UnixTimeToHourMinSec( int rawtime, int &hour, int &min, int &sec, int &day = dummy_int, int &mon = dummy_int, int &year = dummy_int );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Получение текущей даты и времени
/// \return Возвращает текущую дату в строке формата YYYY-MM-DD.HH:mm:ss
///
const std::string CurrentDateTimeToString();

}
}
#endif // SPML_CONVERT_H
/// \}
