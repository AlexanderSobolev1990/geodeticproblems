//----------------------------------------------------------------------------------------------------------------------
///
/// \file       test_geodetic.cpp
/// \brief      Тестирование задач о назначениях
/// \date       20.02.21 - создан
/// \author     Соболев А.А.
///

#ifndef TEST_GEODETIC_H
#define TEST_GEODETIC_H

#include <QtTest>
#include <fstream>

#include <geodesy.h>

using namespace std;

class test_geodetic : public QObject
{
    Q_OBJECT

public:
    test_geodetic();
    ~test_geodetic();

private:
    SPML::Units::TAngleUnit unitAngle = SPML::Units::TAngleUnit::AU_Radian;
    SPML::Units::TRangeUnit unitRange = SPML::Units::TRangeUnit::RU_Meter;

    //
    // Внимание! Правильные ответы на тесты GEOtoRAD и RADtoGEO посчитаны на сфере Sphere6371 и эллипсоиде WGS84 !
    //
    SPML::Geodesy::CEllipsoid sphere_6371 = SPML::Geodesy::Sphere6371;
    SPML::Geodesy::CEllipsoid el_WGS84 = SPML::Geodesy::WGS84;

    static double startPoint[25][2];// начальная точка (lat-lon)
    static double endPoint[25][2];  // конечная точка (lat-lon)
    static double rightAnswerInverseProblem_Sphere6371000[25][3];
    static double rightAnswerDirectProblem_Sphere6371000[25][3];
    static double rightAnswerInverseProblem_WGS84[25][3];
    static double rightAnswerDirectProblem_WGS84[25][3];

private slots:
    void test_GEOtoRAD_RADtoGEO_Sphere_float();
    void test_GEOtoRAD_RADToGEO_Sphere_double();

    void test_GEOtoRAD_RADToGEO_WGS84_float();
    void test_GEOToRAD_RADToGEO_WGS84_double();

    void test_GEOtoRAD_Time_Sphere_float();
    void test_GEOtoRAD_Time_Sphere_double();

    void test_RADtoGEO_Time_Sphere_float();
    void test_RADtoGEO_Time_Sphere_double();

    void test_GEOtoRAD_Time_WGS84_float();
    void test_GEOtoRAD_Time_WGS84_double();

    void test_RADtoGEO_Time_WGS84_float();
    void test_RADtoGEO_Time_WGS84_double();
};


// начальная точка (lat-lon)
double test_geodetic::startPoint[25][2] =
{
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 0.0, 0.0 },
    { 55.7522200, 37.6155600 },
    { 55.7522200, 37.6155600 },
    { 55.7522200, 37.6155600 },
    { 55.7522200, 37.6155600 },
    { 55.7522200, 37.6155600 },
    { 55.7522200, 37.6155600 },
    { 55.7522200, 37.6155600 },
    { 55.7522200, 37.6155600 },
    { 55.7522200, 37.6155600 }
};

// конечная точка (lat-lon)
double test_geodetic::endPoint[25][2] =
{
    { 30.0, 0.0 },
    { 60.0, 0.0 },
    { 30.0, 30.0 },
    { 60.0, 60.0 },
    { 0.0, 30.0 },
    { 0.0, 60.0 },
    { -30.0, 30.0 },
    { -60.0, 60.0 },
    { -30.0, 0.0 },
    { -60.0, 0.0 },
    { -30.0, -30.0 },
    { -60.0, -60.0 },
    { 0.0, -30.0 },
    { 0.0, -60.0 },
    { 30.0, -30.0 },
    { 60.0, -60.0 },
    { 59.9386300, 30.3141300 },     // Saint-Petersburg, Russia
    { 52.5243700, 13.4105300 },     // Berlin, Germany
    { 51.5085300, -0.1257400 },     // London, Great Britain
    { 48.8534100, 2.3488000 },      // Paris, France
    { 41.8919300, 12.5113300 },     // Roma, Italy
    { 64.1354800, -21.8954100 },    // Reykjavik, Iceland
    { 40.7142700, -74.0059700 },    // New-York, USA
    { -34.6131500, -58.3772300 },   // Buenos Aires, Argentina
    { -33.8678500, 151.2073200 }    // Sydney, Australia
};

double test_geodetic::rightAnswerInverseProblem_Sphere6371000[25][3] =
{
    { 3335.847799336762, 0.0, 180.0 },
    { 6671.695598673524, 0.0, 180.0 },
    { 4604.53989281927, 40.893394649131, 229.106605350869 },
    { 8397.717492500104, 26.565051177078, 243.434948822922 },
    { 3335.847799336762, 90.0, 270.0 },
    { 6671.695598673524, 90.0, 270.0 },
    { 4604.53989281927, 139.106605350869, 310.893394649131 },
    { 8397.717492500104, 153.434948822922, 296.565051177078 },
    { 3335.847799336762, 180.0, 0.0 },
    { 6671.695598673524, 180.0, 0.0 },
    { 4604.53989281927, 220.893394649131, 49.106605350869 },
    { 8397.717492500104, 206.565051177078, 63.434948822922 },
    { 3335.847799336762, 270.0, 90.0 },
    { 6671.695598673524, 270.0, 90.0 },
    { 4604.53989281927, 319.106605350869, 130.893394649131 },
    { 8397.717492500104, 333.434948822922, 116.565051177078 },

    { 634.430920264603, 320.181220542521, 133.993227500364 },
    { 1608.171539213493, 267.225012356327, 67.500629141408 },
    { 2500.279541412505, 275.046184204728, 64.249739712143 },
    { 2486.515274556244, 266.942048390901, 58.657801947884 },
    { 2375.139101545627, 240.124400843247, 40.960381828756 },
    { 3305.826756198104, 310.708086770129, 77.933216738675 },
    { 7510.302978483882, 310.317895818262, 34.47939646295 },
    { 13475.866954006504, 253.102223319853, 40.864960821945 },
    { 14496.045468976658, 92.928775544395, 317.398997046878 }
};

double test_geodetic::rightAnswerDirectProblem_Sphere6371000[25][3] =
{
    { 30.0, 0.0, 0.0 },
    { 60.0, 0.0, 0.0 },
    { 30.00000, 30.00000, 49.10661 },
    { 60.00000, 60.00000, 63.43495 },
    { 0.00000, 30.00000, 90.00000 },
    { 0.00000, 60.00000, 90.00000 },
    { -30.00000, 30.00000, 130.89339 },
    { -60.00000, 60.00000, 116.56505 },
    { -30.00000, 0.00000, 180.0 },
    { -60.00000, 0.00000, 180.0 },
    { -30.00000, -30.00000, 229.10661 },
    { -60.00000, -60.00000, 243.43495 },
    { 0.00000, -30.00000, 270.00000 },
    { 0.00000, -60.00000, 270.00000 },
    { 30.00000, -30.00000, 310.89339 },
    { 60.00000, -60.00000, 296.56505 }, // 16

    { 59.93863, 30.31413, 313.99323 },
    { 52.52437, 13.41053, 247.50063 },
    { 51.50853, -0.12574, 244.24974 },
    { 48.85341, 2.34880, 238.65780 }, // 20

    { 41.89193, 12.51133, 220.96038 },
    { 64.13548, -21.89541, 257.93322 },
    { 40.71427, -74.00597, 214.47940 },
    { -34.61315, -58.37723, 220.86496 },
    { -33.86785, 151.20732, 137.39900 } // 25
};

double test_geodetic::rightAnswerInverseProblem_WGS84[25][3] =
{
    { 3320.11338490278, 0.0, 180.0 },
    { 6654.07283255148, 0.0, 180.0 },
    { 4596.22309715941, 41.066728349773, 229.282146152291 },
    { 8389.65383524219, 26.6683364994786, 243.559103061872 },
    { 3339.5847237982, 90.0, 270.0 },
    { 6679.16944759641, 90.0, 270.0 },
    { 4596.22309715941, 138.933271650227, 310.717853847709 },
    { 8389.65383524219, 153.331663500521, 296.440896938128 },
    { 3320.11338490278, 180.0, 0.0 },
    { 6654.07283255148, 180.0, 0.0 },
    { 4596.22309715941, 221.066728349773, 49.2821461522913 },
    { 8389.65383524219, 206.668336499479, 63.5591030618716 },
    { 3339.5847237982, 270.0, 90.0 },
    { 6679.16944759641, 270.0, 90.0 },
    { 4596.22309715941, 318.933271650227, 130.717853847709 },
    { 8389.65383524219, 333.331663500521, 116.440896938128 }, // 16

    { 636.015930287104, 320.126944770617, 133.93894388736 },
    { 1613.33708756292, 267.253437457649, 67.5288169619985 },
    { 2508.31301483815, 275.070800591577, 64.2734524345967 },
    { 2493.92247471802, 266.983302187342, 58.698181706953 }, // 20

    { 2379.36533650612, 240.204880398643, 41.0402378345664 },
    { 3317.37231308908, 310.685756113135, 77.9091383324392 },
    { 7531.172180035, 310.353662800418, 34.492624273887 },
    { 13459.2791448219, 253.278786758906, 40.9713451154643 },
    { 14484.0079887903, 92.7281133803279, 317.323828722842 } //25
};

double test_geodetic::rightAnswerDirectProblem_WGS84[25][3] =
{
    { 30.0, 0.0, 0.0 },
    { 60.0, 0.0, 0.0 },
    { 30.00000, 30.00000, 49.282146445946 },
    { 60.00000, 60.00000, 63.5590965976542 },
    { 0.00000, 30.00000, 90.00000 },
    { 0.00000, 60.00000, 90.00000 },
    { -30.00000, 30.00000, 130.717853554054 },
    { -60.00000, 60.00000, 116.440905302634 },
    { -30.00000, 0.00000, 180.0 },
    { -60.00000, 0.00000, 180.0 },
    { -30.00000, -30.00000, 229.282146445946 },
    { -60.00000, -60.00000, 243.559094697366 },
    { 0.00000, -30.00000, 270.00000 },
    { 0.00000, -60.00000, 270.00000 },
    { 30.00000, -30.00000, 310.717838319556 },
    { 60.00000, -60.00000, 296.440905302634 }, // 16

    { 59.9386300, 30.3141300, 313.938959560895 },
    { 52.5243700, 13.4105300, 247.528820831467 },
    { 51.5085300, -0.1257400, 244.273443839833 },
    { 48.8534100, 2.3488000, 238.698183934423 }, // 20

    { 41.8919300, 12.5113300, 221.040246161089 },
    { 64.1354800, -21.8954100, 257.9091202346 },
    { 40.7142700, -74.0059700, 214.492623699066 },
    { -34.6131500, -58.3772300, 220.971345375481 },
    { -33.8678500, 151.2073200, 137.323827838602 } // 25
};


test_geodetic::test_geodetic()
{
    freopen( "test_geodetic.log", "w", stdout );
}

test_geodetic::~test_geodetic()
{

}


void test_geodetic::test_GEOtoRAD_RADtoGEO_Sphere_float()
{
    std::ofstream out( "test_GeoToRad_RadToGeo_Sphere_float.txt" );

    float latEnd;
    float lonEnd;
    float range;
    float az;
    float az2;
    float eps = 0.001f;

    for( int i = 0; i < 25; i++ ) {
        GEOtoRAD(
            sphere_6371, unitRange, unitAngle,
            static_cast<float>( startPoint[i][0] * SPML::Convert::DgToRdD ),
            static_cast<float>( startPoint[i][1] * SPML::Convert::DgToRdD ),
            static_cast<float>( endPoint[i][0] * SPML::Convert::DgToRdD ),
            static_cast<float>( endPoint[i][1] * SPML::Convert::DgToRdD ),
            range,
            az,
            az2 );
        out << "----------------------------------------------------------------------------------------------" << endl;
        out << "Test #" << ( i + 1 ) << endl;
        out << "StartPoint:" << endl;
        out << startPoint[i][0] << ' ' << startPoint[i][1] << endl;
        out << "EndPoint:" << endl;
        out << endPoint[i][0] << ' ' << endPoint[i][1] << endl;

        out << "Result Inverse Problem:" << endl;
        out << range << ' ' << az << ' ' << az2 << endl;
        QCOMPARE( SPML::Compare::AreEqual( range / 1000.0f, static_cast<float>( rightAnswerInverseProblem_Sphere6371000[i][0] ), eps ), true );
        QCOMPARE( SPML::Compare::AreEqual( az * SPML::Convert::RdToDgF,    static_cast<float>( rightAnswerInverseProblem_Sphere6371000[i][1] ), eps ), true );
        QCOMPARE( SPML::Compare::AreEqual( az2 * SPML::Convert::RdToDgF,   static_cast<float>( SPML::Convert::AngleTo360( rightAnswerInverseProblem_Sphere6371000[i][2] + 180.0, SPML::Units::TAngleUnit::AU_Degree ) ), eps ), true );

        RADtoGEO(
            sphere_6371, unitRange, unitAngle,
            static_cast<float>( startPoint[i][0] * SPML::Convert::DgToRdD ),
            static_cast<float>( startPoint[i][1] * SPML::Convert::DgToRdD ),
            range,
            az,
            latEnd,
            lonEnd,
            az2 );

        out << "Result Direct Problem:" << endl;
        out << latEnd << ' ' << lonEnd << ' ' << az2 << endl;
        QCOMPARE( SPML::Compare::AreEqual( latEnd * SPML::Convert::RdToDgF, static_cast<float>( rightAnswerDirectProblem_Sphere6371000[i][0] ), eps ), true );
        QCOMPARE( SPML::Compare::AreEqual( lonEnd * SPML::Convert::RdToDgF, static_cast<float>( rightAnswerDirectProblem_Sphere6371000[i][1] ), eps ), true );
        //QCOMPARE( SPML::Compare::AreEqual( az2 * SPML::Convert::RdToDgF,    static_cast<float>( To360( rightAnswerDirectProblem_Sphere6371000[i][2] + 180.0, AU_Degree ) ), eps ), true );
        QCOMPARE( SPML::Compare::AreEqual( az2 * SPML::Convert::RdToDgF,    static_cast<float>( rightAnswerDirectProblem_Sphere6371000[i][2] ), eps ), true );
    }
    out.close();
}

void test_geodetic::test_GEOtoRAD_RADToGEO_Sphere_double()
{
    std::ofstream out("test_GeoToRad_RadToGeo_Sphere_double.txt");

    double latEnd;
    double lonEnd;
    double range;
    double az;
    double az2;
    double eps = 0.0001;

    for( int i = 0; i < 25; i++ ) {
        GEOtoRAD(
            sphere_6371, unitRange, unitAngle,
            startPoint[i][0] * SPML::Convert::DgToRdD,
            startPoint[i][1] * SPML::Convert::DgToRdD,
            endPoint[i][0] * SPML::Convert::DgToRdD,
            endPoint[i][1] * SPML::Convert::DgToRdD,
            range,
            az,
            az2 );
        out << "----------------------------------------------------------------------------------------------" << endl;
        out << "Test #" << ( i + 1 ) << endl;
        out << "StartPoint:" << endl;
        out << startPoint[i][0] << ' ' << startPoint[i][1] << endl;
        out << "EndPoint:" << endl;
        out << endPoint[i][0] << ' ' << endPoint[i][1] << endl;

        out << "Result Inverse Problem:" << endl;
        out << range << ' ' << az << ' ' << az2 << endl;
        QCOMPARE( SPML::Compare::AreEqual( range / 1000.0, rightAnswerInverseProblem_Sphere6371000[i][0], eps ), true );
        QCOMPARE( SPML::Compare::AreEqual( az * SPML::Convert::RdToDgD,    rightAnswerInverseProblem_Sphere6371000[i][1], eps ), true );
        QCOMPARE( SPML::Compare::AreEqual( az2 * SPML::Convert::RdToDgD,   SPML::Convert::AngleTo360( rightAnswerInverseProblem_Sphere6371000[i][2] + 180, SPML::Units::TAngleUnit::AU_Degree ), eps ), true );

        RADtoGEO(
            sphere_6371, unitRange, unitAngle,
            startPoint[i][0] * SPML::Convert::DgToRdD,
            startPoint[i][1] * SPML::Convert::DgToRdD,
            range,
            az,
            latEnd,
            lonEnd,
            az2 );

        out << "Result Direct Problem:" << endl;
        out << latEnd << ' ' << lonEnd << ' ' << az2 << endl;
        QCOMPARE( SPML::Compare::AreEqual( latEnd * SPML::Convert::RdToDgD, rightAnswerDirectProblem_Sphere6371000[i][0], eps ), true );
        QCOMPARE( SPML::Compare::AreEqual( lonEnd * SPML::Convert::RdToDgD, rightAnswerDirectProblem_Sphere6371000[i][1], eps ), true );
        QCOMPARE( SPML::Compare::AreEqual( az2 * SPML::Convert::RdToDgD,    rightAnswerDirectProblem_Sphere6371000[i][2], eps ), true );
    }
    out.close();
}

void test_geodetic::test_GEOtoRAD_RADToGEO_WGS84_float()
{
    std::ofstream out( "test_GeoToRad_RadToGeo_WGS84_float.txt" );

    float latEnd;
    float lonEnd;
    float range;
    float az;
    float az2;
    float eps = 0.002f;

    for( int i = 0; i < 25; i++ ) {

        GEOtoRAD(
            el_WGS84, unitRange, unitAngle,
            static_cast<float>( startPoint[i][0] * SPML::Convert::DgToRdD ),
            static_cast<float>( startPoint[i][1] * SPML::Convert::DgToRdD ),
            static_cast<float>( endPoint[i][0] * SPML::Convert::DgToRdD ),
            static_cast<float>( endPoint[i][1] * SPML::Convert::DgToRdD ),
            range,
            az,
            az2 );

        out << "----------------------------------------------------------------------------------------------" << endl;
        out << "Test #" << ( i + 1 ) << endl;
        out << "StartPoint:" << endl;
        out << startPoint[i][0] << ' ' << startPoint[i][1] << endl;
        out << "EndPoint:" << endl;
        out << endPoint[i][0] << ' ' << endPoint[i][1] << endl;
        out << "Result Inverse Problem:" << endl;
        out << range << ' ' << az << ' ' << az2 << endl;
        QCOMPARE( SPML::Compare::AreEqual( range / 1000.0f, static_cast<float>( rightAnswerInverseProblem_WGS84[i][0] ), eps ), true );
        QCOMPARE( SPML::Compare::AreEqual( az * SPML::Convert::RdToDgF,    static_cast<float>( rightAnswerInverseProblem_WGS84[i][1] ), eps ), true );
        QCOMPARE( SPML::Compare::AreEqual( az2 * SPML::Convert::RdToDgF,   static_cast<float>( SPML::Convert::AngleTo360( rightAnswerInverseProblem_WGS84[i][2] + 180, SPML::Units::TAngleUnit::AU_Degree ) ), eps ), true );

        RADtoGEO(
            el_WGS84, unitRange, unitAngle,
            static_cast<float>( startPoint[i][0] * SPML::Convert::DgToRdD ),
            static_cast<float>( startPoint[i][1] * SPML::Convert::DgToRdD ),
            range,
            az,
            latEnd,
            lonEnd,
            az2 );

        out << "Result Direct Problem:" << endl;
        out << latEnd << ' ' << lonEnd << ' ' << az2 << endl;
        QCOMPARE( SPML::Compare::AreEqual( latEnd * SPML::Convert::RdToDgF, static_cast<float>( rightAnswerDirectProblem_WGS84[i][0] ), eps ), true );
        QCOMPARE( SPML::Compare::AreEqual( lonEnd * SPML::Convert::RdToDgF, static_cast<float>( rightAnswerDirectProblem_WGS84[i][1] ), eps ), true );
        QCOMPARE( SPML::Compare::AreEqual( az2 * SPML::Convert::RdToDgF,    static_cast<float>( rightAnswerDirectProblem_WGS84[i][2] ), eps ), true );
    }
    out.close();
}

void test_geodetic::test_GEOToRAD_RADToGEO_WGS84_double()
{
    std::ofstream out("test_GeoToRad_RadToGeo_WGS84_double.txt");

    double latEnd;
    double lonEnd;
    double range;
    double az;
    double az2;
    double eps = 0.002;

    for( int i = 0; i < 25; i++ ) {
        GEOtoRAD(
            el_WGS84, unitRange, unitAngle,
            startPoint[i][0] * SPML::Convert::DgToRdD,
            startPoint[i][1] * SPML::Convert::DgToRdD,
            endPoint[i][0] * SPML::Convert::DgToRdD,
            endPoint[i][1] * SPML::Convert::DgToRdD,
            range,
            az,
            az2 );
        out << "----------------------------------------------------------------------------------------------" << endl;
        out << "Test #" << ( i + 1 ) << endl;
        out << "StartPoint:" << endl;
        out << startPoint[i][0] << ' ' << startPoint[i][1] << endl;
        out << "EndPoint:" << endl;
        out << endPoint[i][0] << ' ' << endPoint[i][1] << endl;

        out << "Result Inverse Problem:" << endl;
        out << range << ' ' << az << ' ' << az2 << endl;
        QCOMPARE( SPML::Compare::AreEqual( range / 1000.0, rightAnswerInverseProblem_WGS84[i][0], eps ), true );
        QCOMPARE( SPML::Compare::AreEqual( az * SPML::Convert::RdToDgD,    rightAnswerInverseProblem_WGS84[i][1], eps ), true );
        QCOMPARE( SPML::Compare::AreEqual( az2 * SPML::Convert::RdToDgD,   SPML::Convert::AngleTo360( rightAnswerInverseProblem_WGS84[i][2] + 180, SPML::Units::TAngleUnit::AU_Degree ), eps ), true );

        RADtoGEO(
            el_WGS84, unitRange, unitAngle,
            startPoint[i][0] * SPML::Convert::DgToRdD,
            startPoint[i][1] * SPML::Convert::DgToRdD,
            range,
            az,
            latEnd,
            lonEnd,
            az2 );

        out << "Result Direct Problem:" << endl;
        out << latEnd << ' ' << lonEnd << ' ' << az2 << endl;
        QCOMPARE( SPML::Compare::AreEqual( latEnd * SPML::Convert::RdToDgD, rightAnswerDirectProblem_WGS84[i][0], eps ), true );
        QCOMPARE( SPML::Compare::AreEqual( lonEnd * SPML::Convert::RdToDgD, rightAnswerDirectProblem_WGS84[i][1], eps ), true );
        QCOMPARE( SPML::Compare::AreEqual( az2 * SPML::Convert::RdToDgD,    rightAnswerDirectProblem_WGS84[i][2], eps ), true );
    }
    out.close();
}

void test_geodetic::test_GEOtoRAD_Time_Sphere_float()
{
    float range;
    float az;
    float azEnd;

    QBENCHMARK {
        GEOtoRAD(
            sphere_6371, unitRange, unitAngle,
            static_cast<float>( startPoint[20][0] * SPML::Convert::DgToRdD ),
            static_cast<float>( startPoint[20][1] * SPML::Convert::DgToRdD ),
            static_cast<float>( endPoint[20][0] * SPML::Convert::DgToRdD ),
            static_cast<float>( endPoint[20][1] * SPML::Convert::DgToRdD ),
            range,
            az,
            azEnd );
    }
}

void test_geodetic::test_GEOtoRAD_Time_Sphere_double()
{
    double range;
    double az;
    double azEnd;

    QBENCHMARK {
        GEOtoRAD(
            sphere_6371, unitRange, unitAngle,
            startPoint[20][0] * SPML::Convert::DgToRdD,
            startPoint[20][1] * SPML::Convert::DgToRdD,
            endPoint[20][0] * SPML::Convert::DgToRdD,
            endPoint[20][1] * SPML::Convert::DgToRdD,
            range,
            az,
            azEnd );
    }
}

void test_geodetic::test_RADtoGEO_Time_Sphere_float()
{
    float latEnd;
    float lonEnd;
    float range;
    float az;
    float azEnd;

    GEOtoRAD(
        sphere_6371, unitRange, unitAngle,
        static_cast<float>( startPoint[20][0] * SPML::Convert::DgToRdD ),
        static_cast<float>( startPoint[20][1] * SPML::Convert::DgToRdD ),
        static_cast<float>( endPoint[20][0] * SPML::Convert::DgToRdD ),
        static_cast<float>( endPoint[20][1] * SPML::Convert::DgToRdD ),
        range,
        az,
        azEnd );

    QBENCHMARK {
        RADtoGEO(
            sphere_6371, unitRange, unitAngle,
            static_cast<float>( startPoint[20][0] * SPML::Convert::DgToRdD ),
            static_cast<float>( startPoint[20][1] * SPML::Convert::DgToRdD ),
            range,
            az,
            latEnd,
            lonEnd,
            azEnd );
    }
}

void test_geodetic::test_RADtoGEO_Time_Sphere_double()
{
    double latEnd;
    double lonEnd;
    double range;
    double az;
    double azEnd;

    GEOtoRAD(
        sphere_6371, unitRange, unitAngle,
        startPoint[20][0] * SPML::Convert::DgToRdD,
        startPoint[20][1] * SPML::Convert::DgToRdD,
        endPoint[20][0] * SPML::Convert::DgToRdD,
        endPoint[20][1] * SPML::Convert::DgToRdD,
        range,
        az,
        azEnd );

    QBENCHMARK {
        RADtoGEO(
            sphere_6371, unitRange, unitAngle,
            startPoint[20][0] * SPML::Convert::DgToRdD,
            startPoint[20][1] * SPML::Convert::DgToRdD,
            range,
            az,
            latEnd,
            lonEnd,
            azEnd );
    }
}

void test_geodetic::test_GEOtoRAD_Time_WGS84_float()
{
    float range;
    float az;
    float azEnd;

    QBENCHMARK {
        GEOtoRAD(
            el_WGS84, unitRange, unitAngle,
            static_cast<float>( startPoint[20][0] * SPML::Convert::DgToRdD ),
            static_cast<float>( startPoint[20][1] * SPML::Convert::DgToRdD ),
            static_cast<float>( endPoint[20][0] * SPML::Convert::DgToRdD ),
            static_cast<float>( endPoint[20][1] * SPML::Convert::DgToRdD ),
            range,
            az,
            azEnd );
    }
}

void test_geodetic::test_GEOtoRAD_Time_WGS84_double()
{
    double range;
    double az;
    double azEnd;

    QBENCHMARK {
        GEOtoRAD(
            el_WGS84, unitRange, unitAngle,
            startPoint[20][0] * SPML::Convert::DgToRdD,
            startPoint[20][1] * SPML::Convert::DgToRdD,
            endPoint[20][0] * SPML::Convert::DgToRdD,
            endPoint[20][1] * SPML::Convert::DgToRdD,
            range,
            az,
            azEnd );
    }
}

void test_geodetic::test_RADtoGEO_Time_WGS84_float()
{
    float latEnd;
    float lonEnd;
    float range;
    float az;
    float azEnd;

    GEOtoRAD(
        el_WGS84, unitRange, unitAngle,
        static_cast<float>( startPoint[20][0] * SPML::Convert::DgToRdD ),
        static_cast<float>( startPoint[20][1] * SPML::Convert::DgToRdD ),
        static_cast<float>( endPoint[20][0] * SPML::Convert::DgToRdD ),
        static_cast<float>( endPoint[20][1] * SPML::Convert::DgToRdD ),
        range,
        az,
        azEnd );

    QBENCHMARK {
        RADtoGEO(
            el_WGS84, unitRange, unitAngle,
            static_cast<float>( startPoint[20][0] * SPML::Convert::DgToRdD ),
            static_cast<float>( startPoint[20][1] * SPML::Convert::DgToRdD ),
            range,
            az,
            latEnd,
            lonEnd,
            azEnd );
    }
}

void test_geodetic::test_RADtoGEO_Time_WGS84_double()
{
    double latEnd;
    double lonEnd;
    double range;
    double az;
    double azEnd;

    GEOtoRAD(
        el_WGS84, unitRange, unitAngle,
        startPoint[20][0] * SPML::Convert::DgToRdD,
        startPoint[20][1] * SPML::Convert::DgToRdD,
        endPoint[20][0] * SPML::Convert::DgToRdD,
        endPoint[20][1] * SPML::Convert::DgToRdD,
        range,
        az,
        azEnd );

    QBENCHMARK {
        RADtoGEO(
            el_WGS84, unitRange, unitAngle,
            startPoint[20][0] * SPML::Convert::DgToRdD,
            startPoint[20][1] * SPML::Convert::DgToRdD,
            range,
            az,
            latEnd,
            lonEnd,
            azEnd );
    }
}

QTEST_APPLESS_MAIN(test_geodetic)

#include "test_geodetic.moc"

#endif
