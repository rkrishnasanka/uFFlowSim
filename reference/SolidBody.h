#ifndef SOLIDBODY_H
#define SOLIDBODY_H

class SolidBody {
protected:
    double _posX;
    double _posY;
    double _scaleX;
    double _scaleY;
    double _theta;
    
    double _velX;
    double _velY;
    double _velTheta;
    
    void globalToLocal(double &x, double &y) const;
    void localToGlobal(double &x, double &y) const;
    
    SolidBody(double posX, double posY, double scaleX, double scaleY,
        double theta, double velX, double velY, double velTheta) :
            _posX(posX), _posY(posY), _scaleX(scaleX), _scaleY(scaleY),
            _theta(theta), _velX(velX), _velY(velY), _velTheta(velTheta) {}
                
    virtual ~SolidBody() {};
    
public:
    virtual double distance(double x, double y) const = 0;
    virtual void closestSurfacePoint(double &x, double &y) const = 0;
    virtual void distanceNormal(double &nx, double &ny, double x, double y) const = 0;
    
    double velocityX(double x, double y) const;
    double velocityY(double x, double y) const;
    void velocity(double &vx, double &vy, double x, double y) const;
    void update(double timestep);
};

class SolidBox: public SolidBody {
public:
    
    SolidBox(double x, double y, double sx, double sy, double t, double vx, double vy, double vt) :
        SolidBody(x, y, sx, sy, t, vx, vy, vt) {}

    double distance(double x, double y) const;
    void closestSurfacePoint(double &x, double &y) const;
    void distanceNormal(double &nx, double &ny, double x, double y) const;
};

class SolidSphere: public SolidBody {
public:
    
    SolidSphere(double x, double y, double sx, double sy, double t, double vx, double vy, double vt) :
        SolidBody(x, y, sx, sy, t, vx, vy, vt) {}
    
    double distance(double x, double y) const;
    void closestSurfacePoint(double &x, double &y) const;
    void distanceNormal(double &nx, double &ny, double x, double y) const;
};

#endif
