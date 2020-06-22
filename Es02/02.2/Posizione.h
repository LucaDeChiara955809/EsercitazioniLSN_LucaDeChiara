#ifndef __Posizione__
#define __Posizione__

class Posizione {

private:

protected:
    double m_x, m_y, m_z;

public:
    // constructors
    Posizione();
    Posizione(double x,double y, double z);
    // destructor
    ~Posizione();
    // methods
    void SetX(double);
    void SetY(double);
    void SetZ(double);
    double GetX();
    double GetY();
    double GetZ();
    double GetDist();
};

#endif // __Posizione__
#pragma once