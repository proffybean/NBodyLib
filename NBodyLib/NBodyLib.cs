/*
 * Copyright  (C) 2011 Jeff Balsley. All Rights Reserved
 * 
 * @File:   This Lib contains the Body class used in my N-body code.  It also contains
 *          the physical and astronomical constants used in the code.  The constants can
 *          be declared in MKS or CGS units. 
 * 
 * @Author: Jeff Balsley
 * @Date: 1/20/2011
 * @Version: 0.9
 * 
 * @TODOs
 * ===============================================================================
 *  - Create Kinetic Energy and Potential Energy properties
 * 
 * @Revision History
 * 
 * Date         Author          Comment
 * ================================================================================
 * 1/20/2011    Jeff            - Created this file
 * 2/4/2011     Jeff            - Added method to return angular momentum
 * 
 */
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Diagnostics;

namespace NBodyLib
{
    public class PhysicalConstants
    {
        //
        // Physical constants
        //
        public readonly double G;                       // Gravitational constant
        public readonly double c;                       // speed of light
        public readonly double e_;                      // Charge of electron

        public readonly double GRAVITATIONAL_CONSTANT;  // Gravitational constant
        public readonly double SPEED_OF_LIGHT;          // speed of light
        public readonly double ELECTRON_CHARGE;         // Charge of electron
        
        public readonly long SECONDS_PER_MINUTE = 60;
        public readonly long SECONDS_PER_HOUR = 60 * 60;
        public readonly long SECONDS_PER_DAY = 24 * 60 * 60;
        public readonly long SECONDS_PER_MONTH = 31536000 / 12;
        public readonly long SECONDS_PER_YEAR = 365 * 24 * 60 * 60;
        public readonly long ONE_MINUTE = 60;
        public readonly long ONE_HOUR = 60 * 60;

        //
        // Astromonical constants
        //
        public readonly double AU;                      // Astronomical Unit
        public readonly double pc;                      // Parsec
        public readonly double ly;                      // light year
        public readonly double Mearth;                  // Mass of Earth
        public readonly double Rearth;                  // Radius of Earth
        public readonly double Msun;                    // Mass of Sun
        public readonly double Mmoon;                   // Mass of moon
        public readonly double Rmoon;                   // Radius of Moon
        public readonly double RmoonOrbit;              // Radius of Moon Orbit

        public readonly double ASTRONOMICAL_UNIT;       // Astronomical Unit
        public readonly double PARSEC;                  // Parsec
        public readonly double LIGHT_YEAR;              // light year
        public readonly double MASS_EARTH;              // Mass of Earth
        public readonly double RADIUS_EARTH;            // Radius of Earth
        public readonly double MASS_SUN;                // Mass of Sun
        public readonly double MASS_MOON;               // Mass of moon
        public readonly double MASS_MERCURY;            // Mass; kg
        public readonly double MASS_VENUS;              // Mass; kg
        public readonly double MASS_MARS;               // Mass; kg
        public readonly double MASS_JUPITER;            // Mass; kg
        public readonly double MERCURY_ORBIT_SPEED;     // Average orbital speed; m/s
        public readonly double VENUS_ORBIT_SPEED;       // Average orbital speed; m/s
        public readonly double EARTH_ORBIT_SPEED;       // Average orbital speed; m/s
        public readonly double MARS_ORBIT_SPEED;        // Average orbital speed; m/s
        public readonly double JUPITER_ORBIT_SPEED;     // Average orbital speed; m/s
        public readonly double MERCURY_ORBIT;           // Ave Orbit; m
        public readonly double VENUS_ORBIT;             // Ave Orbit; m
        public readonly double MARS_ORBIT;              // Ave Orbit; m
        public readonly double JUPITER_ORBIT;           // Ave Orbit; m
        public readonly double RADIUS_MOON;             // Radius of Moon
        public readonly double MOON_ORBIT_RADIUS;       // Radius of Moon Orbit
        public readonly double AVERAGE_MOON_ORBIT_SPEED;// Average orbital speed of Moon

        public readonly double MASS_HALLEYS_COMET;      // Mass of Halley's comet;

        public PhysicalConstants(string units)
        {
            if (String.Equals(units, "CGS", StringComparison.OrdinalIgnoreCase))
            {
                //
                // Physical constants in CGS units
                //
                G = 6.67259 - 8;                    // dyne cm^2 / g^2 Gravitational constant
                c = 2.99792458e10;                  // cm/s speed of light
                e_ = 4.803206e-10;                  // Charge of electron; esu

                GRAVITATIONAL_CONSTANT = 6.67259e-8;// Nm^2 / kg^2 Gravitational constant
                SPEED_OF_LIGHT = 2.99792458e10;     // m/s speed of light
                ELECTRON_CHARGE = 4.803206e-10;     // C Charge of electron

                //
                // Astromonical constants in CGS units
                //
                AU = 1.4960e13;                     // Astronomical Unit; cm
                pc = 3.0857e18;                     // Parsec; cm
                ly = 9.4605e17;                     // light year; cm
                Mearth = 5.974e27;                  // Mass of Earth; g
                Rearth = 6.378e8;                   // Radius of Earth; cm
                Msun = 1.989e33;                    // Mass of Sun; g
                Mmoon = 7.35e25;                    // Mass of moon; g
                Rmoon = 1.738e8;                    // Radius of Moon; cm
                RmoonOrbit = 3.84e10;               // Radius of Moon Orbit; cm

                ASTRONOMICAL_UNIT = 1.4960e13;      // Astronomical Unit; cm
                PARSEC = 3.0857e18;                 // Parsec; cm
                LIGHT_YEAR = 9.4605e17;             // light year; cm
                MASS_EARTH = 5.974e27;              // Mass of Earth; g
                RADIUS_EARTH = 6.378e8;             // Radius of Earth; cm
                MASS_SUN = 1.989e33;                // Mass of Sun; g
                MASS_MOON = 7.35e25;                // Mass of moon; g
                RADIUS_MOON = 1.738e8;              // Radius of Moon; cm
                MOON_ORBIT_RADIUS = 3.84e10;        // Radius of Moon Orbit; cm
                AVERAGE_MOON_ORBIT_SPEED = 102300;  // Average orbital speed of Moon; cm/s
            }
            else if (String.Equals(units, "MKS", StringComparison.OrdinalIgnoreCase))
            {
                G = 6.67e-11; ;                     // Nm^2 / kg^2 Gravitational constant
                c = 2.99792458e8;                   // m/s speed of light
                e_ = 1.60e-19; ;                    // Coulombs  Charge of electron

                GRAVITATIONAL_CONSTANT = 6.67e-11;  // Nm^2 / kg^2 Gravitational constant
                SPEED_OF_LIGHT = 2.99792458e8;      // m/s speed of light
                e_ = 1.60e-19;                      // C Charge of electron

                AU = 1.496e11;                      // Astronomical Unit; m
                pc = 3.086e16;                      // Parsec; m
                ly = 9.46e15;                       // light year; m
                Mearth = 5.98e24;                   // Mass of Earth; kg
                Rearth = 6.378e6;                   // Radius of Earth; m
                Msun = 1.99e30;                     // Mass of Sun; kg
                Mmoon = 7.3477e22;                  // Mass of moon; kg
                Rmoon = 1.738e6;                    // Radius of Moon; m
                RmoonOrbit = 3.84e8;                // Radius of Moon Orbit; m

                ASTRONOMICAL_UNIT = 1.496e11;       // Astronomical Unit; m
                PARSEC = 3.086e16;                  // Parsec; m
                LIGHT_YEAR = 9.46e15;               // light year; m
                MASS_EARTH = 5.98e24;               // Mass of Earth; kg
                RADIUS_EARTH = 6.378e6;             // Radius of Earth; m
                MASS_SUN = 1.99e30;                 // Mass of Sun; kg
                MASS_MOON = 7.3477e22;              // Mass of moon; kg
                MASS_MERCURY = 3.3022e23;           // Mass; kg
                MASS_VENUS = 4.8685e24;             // Mass; kg
                MASS_MARS = 6.4185e23;              // Mass; kg
                MASS_JUPITER = 1.8986e27;           // Mass; kg
                MERCURY_ORBIT_SPEED = 47870;        // Average orbital speed; m/s
                VENUS_ORBIT_SPEED = 35020;          // Average orbital speed; m/s
                EARTH_ORBIT_SPEED = 29780;          // Average orbital speed; m/s
                MARS_ORBIT_SPEED = 24077;           // Average orbital speed; m/s
                JUPITER_ORBIT_SPEED = 13070;        // Average orbital speed; m/s
                MERCURY_ORBIT = 5.79e10;            // Ave Orbit; m
                VENUS_ORBIT = 1.08e11;              // Ave Orbit; m
                MARS_ORBIT = 2.27e11;               // Ave Orbit; m
                JUPITER_ORBIT = 7.78547e11;         // Ave Orbit; m

                RADIUS_MOON = 1.738e6;              // Radius of Moon; m
                MOON_ORBIT_RADIUS = 3.84e8;         // Radius of Moon Orbit; m
                AVERAGE_MOON_ORBIT_SPEED = 1023;    // Average orbital speed of Moon; m/s

                MASS_HALLEYS_COMET = 2.2e14;        // Mass of Halley's comet; kg
            }
            else
            {
                Debug.Assert(false, "Invalid units.  Try passing cgs or mks to constructor");
            }
        
        }
    } // end class PhysicalConstants


    /// <summary>
    /// Author: JBalsley
    /// Date: 1/19/2011 12:14 AM
    /// </summary>
    public class Body
    {
        private string _name;
        public string Name
        {
            set { _name = value; }
            get { return _name; }
        }

        private double _x;
        public double X
        {
            set { _x = value; }
            get { return _x; }
        }

        private double _y;
        public double Y
        {
            set { _y = value; }
            get { return _y; }
        }

        private double _vx;
        public double Vx
        {
            set { _vx = value; }
            get { return _vx; }
        }

        private double _vy;
        public double Vy
        {
            set { _vy = value; }
            get { return _vy; }
        }

        private double _mass;
        public double Mass
        {
            set { _mass = value; }
            get { return _mass; }
        }

        private double _fX;
        public double Fx
        {
            set { _fX = value; }
            get { return _fX; }
        }

        private double _fY;
        public double Fy
        {
            set { _fY = value; }
            get { return _fY; }
        }

        private double _electricCharge;
        public double ElectricCharge
        {
            set { _electricCharge = value; }
            get { return _electricCharge; }
        }

        private double _lastPrintedX;
        public double LastPrintedX
        {
            get { return _lastPrintedX; }
            set { _lastPrintedX = value; }
        }

        private double _lastPrintedY;
        public double LastPrintedY
        {
            get { return _lastPrintedY; }
            set { _lastPrintedY = value; }
        }

        private double _kineticEnergy;
        public double KineticEnergy
        {
            get { return _kineticEnergy; }
            set { _kineticEnergy = value; }
        }

        private double _potentialEnergy;
        public double PotentialEnergy
        {
            get { return _potentialEnergy; }
            set { _potentialEnergy = value; }
        }


        public double TotalVelocity()
        {
            return Math.Sqrt(Math.Pow(_vx, 2) + Math.Pow(_vy, 2));
        }

        public void PrintVelocity()
        {
            Console.Write("Vx= {0:e}\tVy= {1:e}\tV= {2:e}", _vx, _vy, TotalVelocity());
        }

        public void PrintPosition()
        {
            Console.Write("X= {0:e10}\tY= {1:e10}", _x, _y);
        }

        public double DistanceFromSun()
        {
            return Math.Sqrt(Math.Pow(_x, 2) + Math.Pow(_y, 2));
        }

        public double KE()
        {
            _kineticEnergy = .5 * Mass * (Math.Pow(_vx, 2) + Math.Pow(_vy, 2));
            return _kineticEnergy;
        }

        //
        // L = p x r => m(v x r)
        //
        public double AngularMomentum()
        {
            double angle = Math.Abs(Math.Atan2(this.X, this.Y));
            double velocity_perp = this.Vx * Math.Cos(angle) + this.Vy * Math.Sin(angle);
            double L = this.Mass * velocity_perp * Math.Sqrt(Math.Pow(this.X, 2) + Math.Pow(this.Y, 2)); // m * v_perp * r

            return Math.Abs(L);
        }

        // TODO: Create Kinetic Energy and Potential Energy methods / properties

        public Body(double mass,
                    double x,
                    double y,
                    double initial_Vx,
                    double initial_Vy,
                    string name)
        {
            this._mass = mass;
            this._x = x;
            this._y = y;
            this._vx = initial_Vx;
            this._vy = initial_Vy;
            this._name = name;

            // set last printed
            this._lastPrintedX = x;
            this._lastPrintedY = y;
        }

    } // end class Body

    public struct EnergyStruct
    {
        public double Kinetic;
        public double Potential;
        public double Total;

        public EnergyStruct(double KE, double PE, double TE)
        {
            Kinetic = KE;
            Potential = PE;
            Total = TE;
        }
    }
}
