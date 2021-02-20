//----------------------------------------------------------------------------------------------------------------------
///
/// \file       compare.h
/// \brief      Функции сравнения чисел, массивов
/// \date       27.07.20 - создан
/// \author     Соболев А.А.
/// \addtogroup spml
/// \{
///

#ifndef SPML_COMPARE_H
#define SPML_COMPARE_H

// SPML includes:
#include <consts.h>

namespace SPML /// Библиотека СБПМ
{
namespace Compare /// Сравнение чисел
{
//----------------------------------------------------------------------------------------------------------------------
static const float EPS_F = 1.0e-4f; ///< Точность по умолчанию при сравнениях чисел типа float (1.0e-4)
static const float EPS_D = 1.0e-8; ///< Точность по умолчанию при сравнениях чисел типа double (1.0e-8)

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Сравнение двух действительных чисел
/// \param[in] first  - первое число
/// \param[in] second - второе число
/// \param[in] eps - точность сравнения
/// \return true - если разница меньше точности, иначе false
///
inline bool AreEqual( float first, float second, const float &eps = EPS_F )
{
    return ( std::abs( first - second ) < eps );
}

///
/// \brief Сравнение двух действительных чисел
/// \param[in] first  - первое число
/// \param[in] second - второе число
/// \param[in] eps - точность сравнения
/// \return true - если разница меньше точности, иначе false
///
inline bool AreEqual( double first, double second, const double &eps = EPS_D )
{
    return ( std::abs( first - second ) < eps );
}

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Сравнение двух целых массивов по значениям
/// \param[in] first  - первый массив
/// \param[in] second - второй массив
/// \param[in] length - сравниваемая длина массивов
/// \return true - значения совпадают, иначе false
///
inline bool AreEqual( const int *first, const int *second, int length )
{
    for( int i = 0; i < length; i++ ) {
        if( *( first + i ) != *( second + i ) ) {
            return false;
        }
    }
    return true;
}

///
/// \brief Сравнение двух вещественных массивов по значениям
/// \param[in] first - первый массив
/// \param[in] second - второй массив
/// \param[in] length - сравниваемая длина массивов
/// \param[in] eps - точность сравнения вещественных чисел
/// \return true - значения совпадают с заданной точностью, иначе false
///
inline bool AreEqual( const float *first, const float *second, int length, const float &eps = EPS_F )
{
    for( int i = 0; i < length; i++ ) {
        if( !AreEqual( *( first + i ), *( second + i ), eps ) ) {
            return false;
        }
    }
    return true;
}

///
/// \brief Сравнение двух вещественных массивов по значениям
/// \param[in] first - первый массив
/// \param[in] second - второй массив
/// \param[in] length - сравниваемая длина массивов
/// \param[in] eps - точность сравнения вещественных чисел
/// \return true - значения совпадают с заданной точностью, иначе false
///
inline bool AreEqual( const double *first, const double *second, int length, const double &eps = EPS_D )
{
    for( int i = 0; i < length; i++ ) {
        if( !AreEqual( *( first + i ), *( second + i ), eps ) ) {
            return false;
        }
    }
    return true;
}

//----------------------------------------------------------------------------------------------------------------------
///
/// \brief Проверка действительного числа на равенство нулю с заданной точностью
/// \param[in] value - проверяемое число
/// \param[in] eps   - точность
/// \return true - если разница меньше точности, иначе false
///
inline bool IsZero( float value, const float &eps = EPS_F )
{
    return ( std::abs( value ) < eps );
}

///
/// \brief Проверка действительного числа на равенство нулю с заданной точностью
/// \param[in] value - проверяемое число
/// \param[in] eps   - точность
/// \return true - если разница меньше точности, иначе false
///
inline bool IsZero( double value, const double &eps = EPS_D )
{
    return ( std::abs( value ) < eps );
}

}
}
#endif // SPML_COMPARE_H
/// \}

