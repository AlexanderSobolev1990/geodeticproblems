//----------------------------------------------------------------------------------------------------------------------
///
/// \file       geodesy.h
/// \brief      Земные эллипсоиды, геодезические задачи (персчеты координат)
/// \date       06.11.19 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#ifndef SPML_GEODESY_H
#define SPML_GEODESY_H

#include <cassert>

// SPML includes:
#include <compare.h>
#include <convert.h>
#include <units.h>

namespace SPML /// Библиотека СБПМ
{
namespace Geodesy /// Геодезические функции
{
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Земной эллипсоид
///
struct CEllipsoid
{
public:

    ///
    /// \brief Большая полуось эллипсоида (экваториальный радиус)
    /// \return Возвращает большую полуось эллипсоида (экваториальный радиус) в [м]
    ///
    double A() const
    {
        return a;
    }

    ///
    /// \brief Малая полуось эллипсоида (полярный радиус)
    /// \return Возвращает малую полуось эллипсоида (полярный радиус) в [м]
    ///
    double B() const
    {
        return b;
    }

    ///
    /// \brief Cжатие f = ( a - b ) / a
    /// \return Возвращает сжатие
    ///
    double F() const
    {
        return f;
    }

    ///
    /// \brief Обратное сжатие Invf = a / ( a - b )
    /// \return Возвращает обратное сжатие
    ///
    double Invf() const
    {
        return invf;
    }

    ///
    /// \brief Первый эксцентриситет эллипсоида e1 = sqrt( ( a * a ) - ( b * b ) ) / a;
    /// \return Возвращает первый эксцентриситет эллипсоида
    ///
    double EccentricityFirst() const
    {
        return ( sqrt( ( a * a ) - ( b * b ) ) / a );
    }

    ///
    /// \brief Второй эксцентриситет эллипсоида e2 = sqrt( ( a * a ) - ( b * b ) ) / b;
    /// \return Возвращает второй эксцентриситет эллипсоида
    ///
    double EccentricitySecond() const
    {
        return ( sqrt( ( a * a ) - ( b * b ) ) / b );
    }

    ///
    /// \brief Конструктор по умолчанию
    ///
    CEllipsoid();

    ///
    /// \brief Параметрический конструктор эллипсоида
    /// \param[in] semiMajorAxis        - большая полуось (экваториальный радиус)
    /// \param[in] semiMinorAxis        - малая полуось (полярный радиус)
    /// \param[in] inverseFlattening    - обратное сжатие invf = a / ( a - b )
    /// \param[in] isInvfDef            - обратное сжатие задано (малая полуось расчитана из большой и обратного сжатия)
    ///
    CEllipsoid( double semiMajorAxis, double semiMinorAxis, double inverseFlattening, bool isInvfDef );

private: // Доступ к параметрам эллипсоида после его создания не предполагается, поэтому private
    double a;       ///< Большая полуось (экваториальный радиус), [м]
    double b;       ///< Малая полуось (полярный радиус) , [м]
    double invf;    ///< Обратное сжатие invf = a / ( a - b )
    double f;       ///< Сжатие f = ( a - b ) / a
};
//----------------------------------------------------------------------------------------------------------------------
///
///                                      Предопределенные земные эллипсоиды:

///
/// \brief Сфера радиусом 6371000.0 [м] (EPSG::7035)
/// \details Обратное сжатие - бесконечность
///
static const CEllipsoid Sphere6371( 6371000.0, 6371000.0, 0.0, false );

///
/// \brief Сфера радиусом 6378000.0 [м]
/// \details Обратное сжатие - бесконечность
///
static const CEllipsoid Sphere6378( 6378000.0, 6378000.0, 0.0, false );

///
/// \brief Эллипсоид WGS84 (EPSG::7030)
/// \details Главная полуось 6378137.0, обратное сжатие 298.257223563
///
static const CEllipsoid WGS84( 6378137.0, 0.0, 298.257223563, true );

///
/// \brief Эллипсоид GRS80 (EPSG::7019)
/// \details Главная полуось 6378137.0, обратное сжатие 298.257222101
///
static const CEllipsoid GRS80( 6378137.0, 0.0, 298.257222101, true );

///
/// \brief Эллипсоид ПЗ-90 (EPSG::7054)
/// \details Главная полуось 6378136.0, обратное сжатие 298.257839303
///
static const CEllipsoid PZ90( 6378136.0, 0.0, 298.257839303, true );

///
/// \brief Эллипсоид Красовского 1940 (EPSG::7024)
/// \details Главная полуось 6378245.0, обратное сжатие 298.3
///
static const CEllipsoid Krassowsky1940( 6378245.0, 0.0, 298.3, true );

///
/// \brief Сфера радиусом большой полуоси эллипсоида Красовского 1940 (EPSG::7024)
/// \details Обратное сжатие - бесконечность
///
static const CEllipsoid SphereKrassowsky1940( 6378245.0, 6378245.0, 0.0, false );

//----------------------------------------------------------------------------------------------------------------------
///
///                                          Функции пересчета координат
///
static float dummy_float;   // Заглушка для списка параметров функций без перегрузки
static double dummy_double; // Заглушка для списка параметров функций без перегрузки
//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Пересчет географических координат в радиолокационные (Обратная геодезическая задача)
/// \details    Расчет на эллипсоиде по формулам Винсента:
///             \n Vincenty, Thaddeus (April 1975a). "Direct and Inverse Solutions of Geodesics on the Ellipsoid with
///             application of nested equations". Survey Review. XXIII (176): 88–93.
///             \n Расчет на сфере:
///             \n Морозов В.П. Курс сфероидической геодезии. Изд. 2, перераб и доп. М.,Недра, 1979, 296 с., стр 97-100
/// \param[in]  ellipsoid - земной эллипсоид
/// \param[in]  rangeUnit - единицы измерения дальности
/// \param[in]  angleUnit - единицы измерения углов
/// \param[in]  latStart  - широта начальной точки
/// \param[in]  lonStart  - долгота начальной точки
/// \param[in]  latEnd    - широта конечной точки
/// \param[in]  lonEnd    - долгота конечной точки
/// \param[out] d         - расстояние между начальной и конечной точками по ортодроме
/// \param[out] az        - азимут из начальной точки на конечную
/// \param[out] azEnd     - азимут в конечной точке
///
void GEOtoRAD( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    float latStart, float lonStart, float latEnd, float lonEnd, float &d, float &az, float &azEnd = dummy_float );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Пересчет географических координат в радиолокационные (Обратная геодезическая задача)
/// \details    Расчет на эллипсоиде по формулам Винсента:
///             \n Vincenty, Thaddeus (April 1975a). "Direct and Inverse Solutions of Geodesics on the Ellipsoid with
///             application of nested equations". Survey Review. XXIII (176): 88–93.
///             \n Расчет на сфере:
///             \n Морозов В.П. Курс сфероидической геодезии. Изд. 2, перераб и доп. М.,Недра, 1979, 296 с., стр 97-100
/// \param[in]  ellipsoid - земной эллипсоид
/// \param[in]  rangeUnit - единицы измерения дальности
/// \param[in]  angleUnit - единицы измерения углов
/// \param[in]  latStart  - широта начальной точки
/// \param[in]  lonStart  - долгота начальной точки
/// \param[in]  latEnd    - широта конечной точки
/// \param[in]  lonEnd    - долгота конечной точки
/// \param[out] d         - расстояние между начальной и конечной точками по ортодроме
/// \param[out] az        - азимут из начальной точки на конечную
/// \param[out] azEnd     - азимут в конечной точке
///
void GEOtoRAD( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double latStart, double lonStart, double latEnd, double lonEnd, double &d, double &az, double &azEnd = dummy_double );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Пересчет радиолокационных координат в географические (Прямая геодезическая задача)
/// \details    Расчет на эллипсоиде по формулам Винсента:
///             \n Vincenty, Thaddeus (April 1975a). "Direct and Inverse Solutions of Geodesics on the Ellipsoid with
///             application of nested equations". Survey Review. XXIII (176): 88–93.
///             \n Расчет на сфере:
///             \n Морозов В.П. Курс сфероидической геодезии. Изд. 2, перераб и доп. М.,Недра, 1979, 296 с., стр 97-100
/// \param[in]  ellipsoid - земной эллипсоид
/// \param[in]  rangeUnit - единицы измерения дальности
/// \param[in]  angleUnit - единицы измерения углов
/// \param[in]  latStart  - широта начальной точки
/// \param[in]  lonStart  - долгота начальной точки
/// \param[in]  d         - расстояние между начальной и конечной точками по ортодроме
/// \param[in]  az        - азимут из начальной точки на конечную
/// \param[out] latEnd    - широта конечной точки
/// \param[out] lonEnd    - долгота конечной точки
/// \param[out] azEnd     - прямой азимут в конечной точке
///
void RADtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    float latStart, float lonStart, float d, float az, float &latEnd, float &lonEnd, float &azEnd = dummy_float );

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Пересчет радиолокационных координат в географические (Прямая геодезическая задача)
/// \details    Расчет на эллипсоиде по формулам Винсента:
///             \n Vincenty, Thaddeus (April 1975a). "Direct and Inverse Solutions of Geodesics on the Ellipsoid with
///             application of nested equations". Survey Review. XXIII (176): 88–93.
///             \n Расчет на сфере:
///             \n Морозов В.П. Курс сфероидической геодезии. Изд. 2, перераб и доп. М.,Недра, 1979, 296 с., стр 97-100
/// \param[in]  ellipsoid - земной эллипсоид
/// \param[in]  rangeUnit - единицы измерения дальности
/// \param[in]  angleUnit - единицы измерения углов
/// \param[in]  latStart  - широта начальной точки
/// \param[in]  lonStart  - долгота начальной точки
/// \param[in]  d         - расстояние между начальной и конечной точками по ортодроме
/// \param[in]  az        - азимут из начальной точки на конечную
/// \param[out] latEnd    - широта конечной точки
/// \param[out] lonEnd    - долгота конечной точки
/// \param[out] azEnd     - прямой азимут в конечной точке
///
void RADtoGEO( const CEllipsoid &ellipsoid, const Units::TRangeUnit &rangeUnit, const Units::TAngleUnit &angleUnit,
    double latStart, double lonStart, double d, double az, double &latEnd, double &lonEnd, double &azEnd = dummy_double );

}
}
#endif // SPML_GEODESY_H
/// \}
