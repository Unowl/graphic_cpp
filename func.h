
#include <SFML/Graphics.hpp>
#include <string>
#include <vector>
#include <cmath>
#include <limits>
#include <iostream>
#include <iomanip>

template <typename T>
class Point{
private:
    T x_coord;
    T y_coord;

public:
    Point() : x_coord(0), y_coord(0) {}
    Point(T x_, T y_) : x_coord(x_), y_coord(y_) {}
    Point(const Point<T>& p) : x_coord(p.x()), y_coord(p.y()) {}
    T x() const { return x_coord; }
    T y() const { return y_coord; }
    void setx(T x_) { x_coord = x_; }
    void sety(T y_) { y_coord = y_; }
    void print() const { std::cout << x_coord << " " << y_coord << std::endl; }
    Point<T>& operator=(const Point<T>& p) {
        x_coord = p.x();
        y_coord = p.y();
        return *this;
    }
    friend Point<T> operator+(const Point<T>& p1, const Point<T>& p2) { return Point(p1.x() + p2.x(), p1.y() + p2.y()); }
    friend Point<T> operator-(const Point<T>& p1, const Point<T>& p2) { return Point(p1.x() - p2.x(), p1.y() - p2.y()); }
};

template <typename T>
class Point3D {
private:
    T x_coord;
    T y_coord;
    T z_coord;

public:
    Point3D() : x_coord(0), y_coord(0), z_coord(0) {}

    Point3D(T x_, T y_, T z_) : x_coord(x_), y_coord(y_), z_coord(z_) {}

    Point3D(const Point3D<T> &p) : x_coord(p.x()), y_coord(p.y()), z_coord(p.z()) {}

    T x() const { return x_coord; }

    T y() const { return y_coord; }

    T z() const { return z_coord; }

    void setx(T x_) { x_coord = x_; }

    void sety(T y_) { y_coord = y_; }

    void setz(T z_) { z_coord = z_; }

    void print() const { std::cout << x_coord << " " << y_coord << " " << z_coord << std::endl; }

    T LengthSquared() const { return (x_coord * x_coord + y_coord * y_coord + z_coord * z_coord); }

    T Length() const { return (T) sqrt((double) LengthSquared()); }

    Point3D<int> Round() const {
        return Point3D<int>((int) std::round(x_coord), (int) std::round(y_coord), (int) std::round(z_coord));
    }

    Point3D<T> &operator=(const Point3D<T> &p) {
        x_coord = p.x();
        y_coord = p.y();
        z_coord = p.z();
        return *this;
    }

    friend Point3D<T> operator+(const Point3D<T> &p1, const Point3D<T> &p2) {
        return Point3D(p1.x() + p2.x(), p1.y() + p2.y(), p1.z() + p2.z());
    }

    friend Point3D<T> operator-(const Point3D<T> &p1, const Point3D<T> &p2) {
        return Point3D(p1.x() - p2.x(), p1.y() - p2.y(), p1.z() - p2.z());
    }

    const Point3D<T> operator+() const { return Point3D<T>(*this); }

    const Point3D<T> operator-() const { return Point3D<T>(-x_coord, -y_coord, -z_coord); }

    Point3D<T> &operator+=(const Point3D<T> &p) {
        x_coord += p.x();
        y_coord += p.y();
        z_coord += p.z();
        return (*this);
    }

    Point3D<T> &operator-=(const Point3D<T> &p) {
        x_coord -= p.x();
        y_coord -= p.y();
        z_coord -= p.z();
        return (*this);
    }

    Point3D<T> &operator*=(const T &k) {
        x_coord *= k;
        y_coord *= k;
        z_coord *= k;
        return (*this);
    }

    Point3D<T> &operator/=(const T &k) {
        x_coord /= k;
        y_coord /= k;
        z_coord /= k;
        return (*this);
    }

    void Normalize() { *this /= Length(); }

    enum ClassifyResult {
        BEYOND, BETWEEN, BEHIND, RIGHT, LEFT, ORG, DST
    };

    ClassifyResult Classify(const Point3D<T> &org, const Point3D<T> &dst) const;


    void rotate(double phi, Point3D<T> &n) {
        n.Normalize();
        double cos_phi = cos(phi);
        double sin_phi = sin(phi);

        double newx = x_coord * (cos_phi + n.x() * n.x() * (1 - cos_phi)) +
                      y_coord * (n.x() * n.y() * (1 - cos_phi) - n.z() * sin_phi) +
                      z_coord * (n.x() * n.z() * (1 - cos_phi) + n.y() * sin_phi);

        double newy = x_coord * (n.x() * n.y() * (1 - cos_phi) + n.z() * sin_phi) +
                      y_coord * (cos_phi + n.y() * n.y() * (1 - cos_phi)) +
                      z_coord * (n.y() * n.z() * (1 - cos_phi) - n.x() * sin_phi);

        double newz = x_coord * (n.x() * n.z() * (1 - cos_phi) - n.y() * sin_phi) +
                      y_coord * (n.y() * n.z() * (1 - cos_phi) + n.x() * sin_phi) +
                      z_coord * (cos_phi + n.z() * n.z() * (1 - cos_phi));

        x_coord = newx;
        y_coord = newy;
        z_coord = newz;
    }

    void DoIsometric(){
        double newx = (x_coord + z_coord) * 0.707106781186548;
        double newy = (x_coord - z_coord)*0.408248290463863 + y_coord*0.816496580927726;
        double newz = (z_coord - x_coord + y_coord)*0.577350269189626;

        x_coord = newx;
        y_coord = newy;
        z_coord = newz;
    }

    Point3D<T> multiply(double c) const {
        return Point3D<T>(int(round(c * x_coord)), int(round(c * y_coord)), int(round(c * z_coord)));
    };
};

template<typename T>
Point3D<T> cross(const Point3D<T>& a, const Point3D<T>& b) {
    return Point3D<T>(a.y() * b.z() - a.z() * b.y(), -a.x() * b.z() + a.z() * b.x(), a.x() * b.y() - a.y() * b.x());
}

template <> inline
Point3D<int> Point3D<int>::Round() const
{
    return *this;
}

template <class T> inline
const Point3D<T> operator+(const Point3D<T>& l, const Point3D<T>& r)
{
    return Point3D<T>(l.x() + r.x(), l.y() + r.y(), l.z() + r.z());
}

template <class T> inline
const Point3D<T> operator-(const Point3D<T>& l, const Point3D<T>& r)
{
    return Point3D<T>(l.x() - r.x(), l.y() - r.y(), l.z() - r.z());
}

template <class T> inline
T operator*(const Point3D<T>& l, const Point3D<T>& r)
{
    return (l.x() * r.x() + l.y() * r.y() + l.z() * r.z());
}

template <class T> inline
const Point3D<T> operator*(const Point3D<T>& l, const T& r)
{
    return Point3D<T>(l.x() * r, l.y() * r, l.z() * r);
}

template <class T> inline
const Point3D<T> operator*(const T& l, const Point3D<T>& r)
{
    return Point3D<T>(l * r.x(), l * r.y(), l * r.z());
}

template <class T> inline
const Point3D<T> operator/(const Point3D<T>& l, const double& r)
{
    return Point3D<T>(l.x() / r, l.y() / r, l.z() / r);
}

template <class T> inline
const Point3D<T> operator/(const double& l, const Point3D<T>& r)
{
    return Point3D<T>(l / r.x(), l / r.y(), l / r.z());
}

template <class T> inline
bool operator==(const Point3D<T>& l, const Point3D<T>& r)
{
    return (l.x() == r.x() && l.y() == r.y() && l.z() == r.z());
}

template <class T> inline
bool operator!=(const Point3D<T>& l, const Point3D<T>& r)
{
    return !(l == r);
}

template <class T>
typename Point3D<T>::ClassifyResult Point3D<T>::Classify(const Point3D<T>& org, const Point3D<T>& dst) const
{
    const T epsilon = (T)1 / 1000000;

    const Point3D<T> a = dst - org,
            b = (*this) - org;
    T s = a.x * (-b.y) + a.y * b.x;
    if (Abs(s) < epsilon)    // kludge
    {
        if (a.x * b.x < (T)0 || a.y * b.y < (T)0)
            return BEHIND;
        T la = a.LengthSquared(),
                lb = b.LengthSquared();
        if (lb > la)
            return BEYOND;
        if (Abs(lb) < epsilon)
            return ORG;
        if (Abs(lb - la) < epsilon)
            return DST;
        return BETWEEN;
    }
    if (s > (T)0)
        return LEFT;
    else // s < 0
        return RIGHT;
}

struct Front {
    std::vector<Point3D<double>> points;
    Point3D<double> center;
    Point3D<double> n;

    Front() = default;

    Front(const std::vector<Point3D<double>>& points,
          const Point3D<double>& center,
          const Point3D<double>& n) : points(points), center(center), n(n) {}
};

class Cube {

public:
    Point3D<double> center;
    std::vector<Front> arr;

    Cube(const Point3D<double>& p_min, int a, int b, int h) {
        arr.resize(6);
        std::vector<Point3D<double>> low(4), high(4);
        low[0] = p_min;
        low[1] = { p_min.x(), p_min.y() + b, p_min.z() };
        low[2] = { p_min.x() + a, p_min.y() + b, p_min.z() };
        low[3] = { p_min.x() + a, p_min.y(), p_min.z() };

        for (size_t i = 0; i < low.size(); i++)
            high[i] = { low[i].x(), low[i].y(), low[i].z() + h };

        arr[0] = Front(low, {}, { 0, 0, 1 });
        arr[1] = Front(high, {}, { 0, 0, -1 });

        std::vector<Point3D<double>> left = { low[0], low[1], high[1], high[0] };
        std::vector<Point3D<double>> up = { low[1], low[2], high[2], high[1] };
        std::vector<Point3D<double>> right = { low[2], low[3], high[3], high[2] };
        std::vector<Point3D<double>> down = { low[3], low[0], high[0], high[3] };

        arr[2] = Front(left, {}, { 1, 0, 0 });
        arr[3] = Front(up, {}, { 0, -1, 0 });
        arr[4] = Front(right, {}, { -1, 0, 0 });
        arr[5] = Front(down, {}, { 0, 1, 0 });

        center = Point3D<double>( p_min.x() + a / 2, p_min.y() + b / 2, p_min.z() + h / 2 );
        for (auto& face : arr) {
            Point3D<double> face_center(0, 0, 0);
            for (auto& v : face.points)
                face_center += v;
            face.center =  face_center / 4.0;
        }

        FixNormals();
    }

    Point3D<double> GetCenter() const {
        return center;
    }

    void rotate(double phi, Point3D<double> n) {
        for (auto& face : arr) {
            for (auto& v : face.points) {
                v.rotate(phi, n);
            }
            face.center.rotate(phi, n);
        }
        center.rotate(phi, n);

        FixNormals();
    }

    void DoIsometric() {
        for (auto& face : arr) {
            for (auto& v : face.points) {
                v.DoIsometric();
            }
            face.center.DoIsometric();
        }
        center.DoIsometric();

        FixNormals();
    }

    void FixNormals() {
        for (auto& face : arr) {
            face.n = center - face.center;
        }
    }

    Point3D<double> OnePointTransform(const Point3D<double>& point, double r) const {
        return point.multiply(1.0 / (1 + r * point.z()));
    };

};


template <typename T>
int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

enum class pointType { LEFT, RIGHT, BEHIND, BETWEEN, ORIGIN, DESTINATION };
enum class interceptionType { SAME, PARALLEL, CROSS, NO_CROSS };
enum class pointTypeToPolygonEdge { TOUCHING, CROSS_LEFT, CROSS_RIGHT, INESSENTIAL };
enum class pointToPolygonType { INSIDE, OUTSIDE };
enum class methods { EO, NZW };
enum class clockWiseType { CW, CCW, NONCONVEX };

class WinInstance {
private:
    sf::RenderWindow& window;
    sf::Texture texture;
    sf::Sprite sprite;
    sf::Image image;

    template <typename T>
    interceptionType linesInterception(const Point<T>& a1, const Point<T>& a2, const Point<T>& b1, const Point<T>& b2, double& t);

    pointToPolygonType fillPointEvenOdd(const Point<int>& p, const std::vector<Point<int>>& polygon);

    pointToPolygonType fillPointNonZeroWinding(const Point<int>& p, const std::vector<Point<int>>& polygon);

public:
    WinInstance(sf::RenderWindow&);
    void display(){
        window.display();
    }
    bool isOpen() const{
        return window.isOpen();
    }
    bool pollEvent(sf::Event& event){
        return window.pollEvent(event);
    }
    void close(){
        window.close();
    }
    void drawImage(bool update = true){
        if (update) {
            window.clear(sf::Color::White);
            texture.update(image);
            sprite.setTexture(texture);
            window.draw(sprite);
            //window.clear(sf::Color::Red);
            window.display();
            image.create(window.getSize().x, window.getSize().y, sf::Color::White);
        } else {
            window.clear(sf::Color::White);
            window.draw(sprite);
            window.display();
        }
    }
    void clear(sf::Color color){
        window.clear(color);
    }

    void lineBresenham(const Point<int>& s, const Point<int>& e, const sf::Color& color);

    void lineBresenham(const Point3D<double>& s, const Point3D<double>& e, const sf::Color& color){
        Point<int> s_ = Point((int)std::round(s.x()), (int)std::round(s.y()));
        Point<int> e_ = Point((int)std::round(e.x()), (int)std::round(e.y()));
        lineBresenham(s_, e_, color);

    }

    void HatchedLine(const Point<int>& s, const Point<int>& e, const sf::Color& color);

    void polygon(const std::vector<Point<int>>& vertex, const sf::Color& color);

    template <typename T>
    pointType pointPositionToLineSegment(const Point<T>& p, const Point<T>& start, const Point<T>& end);

    bool checkConvex(const std::vector<Point<int>>& vertex);

    template <typename T>
    interceptionType linesInterceptionCoords(const Point<T>& a1, const Point<T>& a2, const Point<T>& b1, const Point<T>& b2, Point<double>& crossCoords);

    template <typename T>
    interceptionType lineSegmentsInterception(const Point<T>& a1, const Point<T>& a2, const Point<T>& b1, const Point<T>& b2, Point<double>& crossCoords);

    bool checkPolygonSimplicity(const std::vector<Point<int>>& vertex, Point<double>& p);

    void boundingBox(const std::vector<Point<int>>& polygon, std::vector<Point<int>>& res);

    pointTypeToPolygonEdge pttpe(const Point<int>& p, const Point<int>& p1, const Point<int>& p2);

    void fillPolygon(const std::vector<Point<int>>& polygon, const methods& method, const sf::Color& color);

    bool saveImage(const std::string& filename){
        return image.saveToFile(filename);
    }

    void curveBezier3(const std::vector<Point<int>>& points, const sf::Color& color);

    void curveBezier5(const std::vector<Point<int>>& points, const sf::Color& color);

    clockWiseType checkClockWise(const std::vector<Point<int>>& vertex);

    bool clipLineCyrusBeck(const std::vector<Point<int>>& polygon, const Point<int>& p1, const Point<int>& p2, Point<int>& p1_new, Point<int>& p2_new);

    void CatmullRom(std::vector<Point<int>>& points, const sf::Color& color);
    void CatmullRomLines(std::vector<Point<int>>& points, const sf::Color& color);

    void ChangeColorsLine(int x_left, int x_right, int y, const sf::Color& new_color);
    void ChangeColorsPixel(int x, int y, const sf::Color& old_color, const sf::Color& new_color, int N);
    void Zatravka(const sf::Color &old_color, const sf::Color &new_color);

    void DrawBounds(Front& front, const sf::Color& color);
    void DrawCube(Cube& cube, const sf::Color& color);
    void DrawBounds(Cube& cube,const sf::Color& color);
    void DrawOnePointProjection(Cube& cube, double r, const sf::Color& color);
};

WinInstance::WinInstance(sf::RenderWindow& window_) : window(window_) {
    window.clear(sf::Color::White);
    image.create(window.getSize().x, window.getSize().y, sf::Color::White);
    texture.loadFromImage(image);
    sprite.setTexture(texture);
    sprite.setPosition(0, 0);
}


template <typename T>
pointType WinInstance::pointPositionToLineSegment(const Point<T>& p, const Point<T>& start, const Point<T>& end) {
    T ax = end.x() - start.x();
    T by = p.y() - start.y();
    T bx = p.x() - start.x();
    T ay = end.y() - start.y();
    double s = ax * by - bx * ay;
    if (s > 0)
        return pointType::RIGHT;
    if (s < 0)
        return pointType::LEFT;
    if ((ax * bx < 0) || (ay * by < 0))
        return pointType::BEHIND;
    if ((ax * ax + ay * ay) < (bx * bx + by * by))
        return pointType::BEHIND;
    if (start.x() == p.x() && start.y() == p.y())
        return pointType::ORIGIN;
    if (end.x() == p.x() && end.y() == p.y())
        return pointType::DESTINATION;
    return pointType::BETWEEN;
}

template <typename T>
interceptionType WinInstance::linesInterception(const Point<T>& a1, const Point<T>& a2, const Point<T>& b1, const Point<T>& b2, double& t) {
    double nx = b2.y() - b1.y();
    double ny = b1.x() - b2.x();
    double denom = nx * (a2.x() - a1.x()) + ny * (a2.y() - a1.y());
    if (denom == 0) {
        pointType type = pointPositionToLineSegment(b1, b2, a1);
        if (type == pointType::LEFT || type == pointType::RIGHT)
            return interceptionType::PARALLEL;
        else
            return interceptionType::SAME;
    }
    double num = nx * (a1.x() - b1.x()) + ny * (a1.y() - b1.y());
    t = -num / denom;
    return interceptionType::CROSS;
}

template <typename T>
interceptionType WinInstance::linesInterceptionCoords(const Point<T>& a1, const Point<T>& a2, const Point<T>& b1, const Point<T>& b2, Point<double>& crossCoords) {
    double t(0.0);
    interceptionType type = linesInterception(a1, a2, b1, b2, t);
    if (type == interceptionType::PARALLEL || type == interceptionType::SAME) {
        return type;
    }
    crossCoords.setx(a1.x() + t * (a2.x() - a1.x()));
    crossCoords.sety(a1.y() + t * (a2.y() - a1.y()));
}

template <typename T>
interceptionType WinInstance::lineSegmentsInterception(const Point<T>& a1, const Point<T>& a2, const Point<T>& b1, const Point<T>& b2, Point<double>& crossCoords) {
    double t(0.0);
    interceptionType type = linesInterception(b1, b2, a1, a2, t);
    if (type == interceptionType::SAME || type == interceptionType::PARALLEL)
        return type;
    if (t < 0.0 || t > 1.0)
        return interceptionType::NO_CROSS;
    linesInterception(a1, a2, b1, b2, t);
    if (t < 0.0 || t > 1.0)
        return interceptionType::NO_CROSS;
    crossCoords.setx(a1.x() + t * (a2.x() - a1.x()));
    crossCoords.sety(a1.y() + t * (a2.y() - a1.y()));
    return interceptionType::CROSS;
}

void WinInstance::lineBresenham(const Point<int>& s, const Point<int>& e, const sf::Color& color) {
    Point<int> start(s.x(), s.y());
    Point<int> end(e.x(), e.y());
    if (start.x() > end.x())
        std::swap(start, end);
    int x, y;
    int ix, iy;
    int dx = abs(end.x() - start.x()), dy = abs(end.y() - start.y());
    bool swapped(false);

    if (dx >= dy) {
        x = start.x(), y = start.y();
        ix = sgn(end.x() - start.x()), iy = sgn(end.y() - start.y());
    } else {
        x = start.y(), y = start.x();
        ix = sgn(end.y() - start.y()), iy = sgn(end.x() - start.x());
        std::swap(dx, dy);
        swapped = true;
    }

    int error = 2 * dy - dx;
    for (int i = 0; i <= dx; i++) {
        if (!swapped)
            image.setPixel(x, y, color);
        else
            image.setPixel(y, x, color);
        if (error >= 0) {
            y += iy;
            error -= 2 * dx;
        }
        x += ix;
        error += 2 * dy;
    }
}

void WinInstance::polygon(const std::vector<Point<int>>& vertex, const sf::Color& color) {
    for (int i = 0; i < vertex.size() - 1; ++i) {
        lineBresenham(vertex[i], vertex[i + 1], color);
    }
}

bool WinInstance::checkConvex(const std::vector<Point<int>>& vertex) {
    pointType type = pointPositionToLineSegment(vertex[0], vertex[1], vertex[2]);
    for (int i = 0; i < vertex.size() - 2; i++) {
        for (int j = 0; j < vertex.size() - 1; j++) {
            if (j != i && j != i + 1) {
                if (pointPositionToLineSegment(vertex[j], vertex[i], vertex[i + 1]) != type)
                    return false;
            }
        }
    }
    return true;
}

clockWiseType WinInstance::checkClockWise(const std::vector<Point<int>>& vertex) {
    if (!checkConvex(vertex))
        return clockWiseType::NONCONVEX;
    pointType type = pointPositionToLineSegment(vertex[0], vertex[1], vertex[2]);
    if (type == pointType::RIGHT)
        return clockWiseType::CW;
    else
        return clockWiseType::CCW;
}

bool WinInstance::checkPolygonSimplicity(const std::vector<Point<int>>& vertex, Point<double>& p) {
    for (int i = 0; i < vertex.size() - 2; i++) {
        for (int j = 0; j < vertex.size() - 1; j++) {
            if (j != i && j != i + 1 && (j != (i > 0 ? i - 1 : vertex.size() - 2))) {
                if (lineSegmentsInterception(vertex[i], vertex[i + 1], vertex[j], vertex[j + 1], p) == interceptionType::CROSS)
                    return false;
            }
        }
    }
    return true;
}

void WinInstance::boundingBox(const std::vector<Point<int>>& polygon, std::vector<Point<int>>& res) {
    int maxx = 0.0, maxy = 0.0, minx = std::numeric_limits<int>::max(), miny = std::numeric_limits<int>::max();
    for (Point<int> vertex : polygon) {
        if (vertex.x() > maxx)
            maxx = vertex.x();
        if (vertex.x() < minx)
            minx = vertex.x();
        if (vertex.y() > maxy)
            maxy = vertex.y();
        if (vertex.y() < miny)
            miny = vertex.y();
    }
    res[0].setx(minx), res[0].sety(miny);
    res[1].setx(minx), res[1].sety(maxy);
    res[2].setx(maxx), res[2].sety(maxy);
    res[3].setx(maxx), res[3].sety(miny);
    res[4].setx(minx), res[4].sety(miny);
    res.resize(5);
    res.shrink_to_fit();
}

pointTypeToPolygonEdge WinInstance::pttpe(const Point<int>& p, const Point<int>& p1, const Point<int>& p2) {
    switch (pointPositionToLineSegment(p, p1, p2)) {
        case pointType::LEFT:
            if (p.y() > p1.y() && p.y() <= p2.y())
                return pointTypeToPolygonEdge::CROSS_LEFT;
            else
                return pointTypeToPolygonEdge::INESSENTIAL;
        case pointType::RIGHT:
            if (p.y() > p2.y() && p.y() <= p1.y())
                return pointTypeToPolygonEdge::CROSS_RIGHT;
            else
                return pointTypeToPolygonEdge::INESSENTIAL;
        case pointType::BETWEEN:
        case pointType::ORIGIN:
        case pointType::DESTINATION:
            return pointTypeToPolygonEdge::TOUCHING;
        default:
            return pointTypeToPolygonEdge::INESSENTIAL;
    }
}

pointToPolygonType WinInstance::fillPointEvenOdd(const Point<int>& p, const std::vector<Point<int>>& polygon) {
    int param(0);
    for (int i = 0; i < polygon.size() - 1; i++) {
        switch (pttpe(p, polygon[i], polygon[i + 1])) {
            case pointTypeToPolygonEdge::TOUCHING:
                return pointToPolygonType::INSIDE;
            case pointTypeToPolygonEdge::CROSS_LEFT:
            case pointTypeToPolygonEdge::CROSS_RIGHT:
                param = 1 - param;
                break;
        }
    }
    if (param)
        return pointToPolygonType::INSIDE;
    else
        return pointToPolygonType::OUTSIDE;
}

pointToPolygonType WinInstance::fillPointNonZeroWinding(const Point<int>& p, const std::vector<Point<int>>& polygon) {
    int param(0);
    for (int i = 0; i < polygon.size() - 1; i++) {
        switch (pttpe(p, polygon[i], polygon[i + 1])) {
            case pointTypeToPolygonEdge::TOUCHING:
                return pointToPolygonType::INSIDE;
            case pointTypeToPolygonEdge::CROSS_LEFT:
                param++;
                break;
            case pointTypeToPolygonEdge::CROSS_RIGHT:
                param--;
                break;
        }
    }
    if (param)
        return pointToPolygonType::INSIDE;
    else
        return pointToPolygonType::OUTSIDE;
}

void WinInstance::fillPolygon(const std::vector<Point<int>>& polygon, const methods& method, const sf::Color& color) {
    std::vector<Point<int>> boundbox(5);
    boundingBox(polygon, boundbox);
    switch (method) {
        case methods::EO:
            for (int i = boundbox[0].x(); i <= boundbox[2].x(); i++) {
                for (int j = boundbox[0].y(); j <= boundbox[2].y(); j++) {
                    if (fillPointEvenOdd(Point(i, j), polygon) == pointToPolygonType::INSIDE)
                        image.setPixel(i, j, color);
                }
            }
            break;
        case methods::NZW:
            for (int i = boundbox[0].x(); i <= boundbox[2].x(); i++) {
                for (int j = boundbox[0].y(); j <= boundbox[2].y(); j++) {
                    if (fillPointNonZeroWinding(Point(i, j), polygon) == pointToPolygonType::INSIDE)
                        image.setPixel(i, j, color);
                }
            }
            break;
    }
}


void WinInstance::curveBezier3(const std::vector<Point<int>>& points, const sf::Color& color) {
    if (points.size() != 4) {
        throw std::runtime_error("\n(curveBezier3) Cubic Bézier curves are plotted for 4 points.");
    }
    int H = std::max(abs(points[0].x() - 2 * points[1].x() + points[2].x()) + abs(points[0].y() - 2 * points[1].y() + points[2].y()),
                     abs(points[1].x() - 2 * points[2].x() + points[3].x()) + abs(points[1].y() - 2 * points[2].y() + points[3].y()));
    int N = std::ceil(1 + std::sqrt(3 * H));
    double tau = 1.0 / N;
    double t = tau;
    Point<int> p1, p2(points[0]);
    for (int i = 0; i < N - 1; ++i) {
        p1 = p2;
        p2.setx(std::round(std::pow(1 - t, 3) * points[0].x() + 3 * t * std::pow(1 - t, 2) * points[1].x() + 3 * std::pow(t, 2) * (1 - t) * points[2].x() +
                           std::pow(t, 3) * points[3].x()));
        p2.sety(std::round(std::pow(1 - t, 3) * points[0].y() + 3 * t * std::pow(1 - t, 2) * points[1].y() + 3 * std::pow(t, 2) * (1 - t) * points[2].y() +
                           std::pow(t, 3) * points[3].y()));
        t += tau;
        lineBresenham(p1, p2, color);
    }
    lineBresenham(p2, points[3], color);
}

void WinInstance::curveBezier5(const std::vector<Point<int>>& points, const sf::Color& color) {
    if (points.size() != 6) {
        throw std::runtime_error("\n(curveBezier3) Cubic Bézier curves are plotted for 6 points.");
    }
    //int H = std::max(abs(points[0].x() - 2 * points[1].x() + points[2].x()) + abs(points[0].y() - 2 * points[1].y() + points[2].y()),
    //                 abs(points[1].x() - 2 * points[2].x() + points[3].x()) + abs(points[1].y() - 2 * points[2].y() + points[3].y()));
    //int N = std::ceil(1 + std::sqrt(3 * H));
    int N = 50;
    double tau = 1.0 / N;
    double t = tau;
    Point<int> p1, p2(points[0]);
    for (int i = 0; i < N - 1; ++i) {
        p1 = p2;
        p2.setx(std::round(std::pow(1 - t, 5) * points[0].x() + 5 * t * std::pow(1 - t, 4) * points[1].x() + 10 * std::pow(t, 2) * std::pow(1 - t, 3)  * points[2].x() +
                10 * std::pow(t, 3) * std::pow(1 - t, 2) * points[3].x() + 5 * (1-t) * std::pow(t, 4) * points[4].x() + std::pow(t, 5) * points[5].x()));
        p2.sety(std::round(std::pow(1 - t, 5) * points[0].y() + 5 * t * std::pow(1 - t, 4) * points[1].y() + 10 * std::pow(t, 2) * std::pow(1 - t, 3)  * points[2].y() +
                           10 * std::pow(t, 3) * std::pow(1 - t, 2) * points[3].y() + 5 * (1-t) * std::pow(t, 4) * points[4].x() + std::pow(t, 5) * points[5].y()));
        t += tau;
        lineBresenham(p1, p2, color);
    }
    lineBresenham(p2, points[5], color);
}

bool WinInstance::clipLineCyrusBeck(const std::vector<Point<int>>& polygon, const Point<int>& p1, const Point<int>& p2, Point<int>& p1_new, Point<int>& p2_new) {
    Point<double> __tmp__;
    std::vector<Point<int>> _polygon(polygon.begin(), polygon.end());
    if (!checkPolygonSimplicity(polygon, __tmp__)) {
        throw std::runtime_error("\n(clipLineCyrusBeck) Polygon is not simple.");
    }
    switch (checkClockWise(polygon)) {
        case clockWiseType::NONCONVEX:
            throw std::runtime_error("\n(clipLineCyrusBeck) Polygon is not convex.");
        case clockWiseType::CCW:
            std::reverse(_polygon.begin(), _polygon.end());
            break;
        case clockWiseType::CW:
            break;
    }

    Point<int> s = p2 - p1;
    double t(0.0), t1(0.0), t2(1.0);
    for (int i = 0; i < _polygon.size() - 1; i++) {
        switch (linesInterception(p1, p2, _polygon[i], _polygon[i + 1], t)) {
            case interceptionType::SAME:
                return false;
            case interceptionType::PARALLEL:
                if (pointPositionToLineSegment(p1, _polygon[i], _polygon[i + 1]) == pointType::LEFT)
                    return false;
                break;
            case interceptionType::CROSS:
                int nx = _polygon[i].y() - _polygon[i + 1].y();
                int ny = _polygon[i + 1].x() - _polygon[i].x();
                if (nx * s.x() + ny * s.y() > 0) {
                    t1 = std::max(t, t1);
                } else {
                    t2 = std::min(t, t2);
                }
                break;
        }
    }
    if (t1 <= t2) {
        p1_new.setx(std::round(p1.x() + t1 * (p2.x() - p1.x())));
        p1_new.sety(std::round(p1.y() + t1 * (p2.y() - p1.y())));
        p2_new.setx(std::round(p1.x() + t2 * (p2.x() - p1.x())));
        p2_new.sety(std::round(p1.y() + t2 * (p2.y() - p1.y())));
        return true;
    }
    return false;
}

void WinInstance::HatchedLine(const Point<int>& s, const Point<int>& e, const sf::Color& color) {

    int l_h_x = abs(s.x() - e.x()) / 10;
    int l_h_y = abs(s.y() - e.y()) / 10;
    int new_1_x, new_1_y;
    int new_2_x = s.x();
    int new_2_y = s.y();

    for (int i = 0; i<5; i++){
        new_1_x = new_2_x;
        new_1_y = new_2_y;
        new_2_x = new_1_x + l_h_x;
        new_2_y = new_1_y + l_h_y;

        lineBresenham(Point(new_1_x, new_1_y), Point(new_2_x, new_2_y), color);

        new_1_x = new_2_x;
        new_1_y = new_2_y;
        new_2_x = new_1_x + l_h_x;
        new_2_y = new_1_y + l_h_y;


    }

}

void WinInstance::CatmullRom(std::vector<Point<int>>& points, const sf::Color& color) {
    int dist_1 = abs(points[0].x() - 2 * points[1].x() + points[2].x()) + abs(points[0].y() - 2 * points[1].y() + points[2].y());
    int dist_2 = abs(points[1].x() - 2 * points[2].x() + points[3].x()) + abs(points[1].y() - 2 * points[2].y() + points[3].y());
    int H = (dist_1 > dist_2) ? dist_1 : dist_2;
    int N = std::round(1 + std::sqrt(3 * H));
    double t0 = 1.0 / N;
    double t = 0;
    Point<int>  p1;
    Point<int> p2 = points[1];
    for (int i = 0; i < N - 1; i++) {
        p1 = p2;
        t = t0 * (i + 1);
        p2.setx(std::round(1.0 / 2 * (std::pow(1 - t, 2) * (-t) * points[0].x() + (2 - 5 * t * t + 3 * std::pow(t, 3)) * points[1].x() + t * (1 + 4 * t - 3 * t * t) * points[2].x() -
                                      std::pow(t, 2) * (1 - t) * points[3].x())));
        p2.sety(std::round(1.0 / 2 * (std::pow(1 - t, 2) * (-t) * points[0].y() + (2 - 5 * t * t + 3 * std::pow(t, 3)) * points[1].y() + t * (1 + 4 * t - 3 * t * t) * points[2].y() -
                                      std::pow(t, 2) * (1 - t) * points[3].y())));
        lineBresenham(p1, p2, color);
    }
    lineBresenham(p2, points[2], color);
}

void WinInstance::CatmullRomLines(std::vector<Point<int>>& points, const sf::Color& color) {
    int n = points.size();
    if (n < 4) {
        throw std::runtime_error("\n(CatmullRomLines) Задано меньше 4 точек");
        return;
    }
    int count = 4;
    while (count <= n) {
        std::vector<Point<int>> p = { points[count - 4], points[count - 3], points[count - 2], points[count - 1] };
        CatmullRom(p, color);
        count++;
    }
}

void WinInstance::ChangeColorsLine(int x_left, int x_right, int y, const sf::Color& new_color) {
    for (int i = x_left; i <= x_right; i++) {
        image.setPixel(i, y, new_color);
    }
}

void WinInstance::ChangeColorsPixel(int x, int y, const sf::Color& old_color, const sf::Color& new_color, int N) {
    std::vector<Point<int>> q;
    q.push_back(Point<int>(x, y));
    while (!q.empty())
    {

        Point<int> p = q.back();
        q.pop_back();
        x = p.x();
        y = p.y();

        int x_left = x, x_right = x;
        while (x_left > 0) {
            x_left--;
            if (image.getPixel(x_left, y) != old_color) {
                x_left++;
                break;
            }
        }

        while (x_right < (N - 1)) {
            x_right++;
            if (image.getPixel(x_right, y) != old_color) {
                x_right--;
                break;
            }
        }

        ChangeColorsLine(x_left, x_right, y, new_color);

        int start = 1;
        int x_min = x_left > 0 ? x_left - 1 : x_left;
        int x_max = x_right < N - 2 ? x_right + 1 : x_right;

        if (y > 0) {
            for (int i = x_min; i <= x_max; i++) {
                if (image.getPixel(i, y - 1) == old_color) {
                    if (start) {
                        q.push_back(Point<int>(i, y - 1));
                        start = 0;
                    }
                }
                else {
                    start = 1;
                }
            }
        }
        start = 1;
        if (y < N - 1) {
            for (int i = x_min; i <= x_max; i++) {
                if (image.getPixel(i, y + 1) == old_color) {
                    if (start) {
                        q.push_back(Point<int>(i, y + 1));
                        start = 0;
                    }
                }
                else {
                    start = 1;
                }
            }
        }
    }
}

void WinInstance::Zatravka(const sf::Color& old_color, const sf::Color& new_color) {
    int N = image.getSize().y;
    for (int x = 0; x < N; x++) {
        for (int y = 0; y < N; y++) {
            if (image.getPixel(x, y) == old_color) {
                ChangeColorsPixel(x, y, old_color, new_color, N);
                return;
            }
        }
    }
}

void WinInstance::DrawBounds(Front& front, const sf::Color& color) {
    lineBresenham(front.points[0], front.points[1], color);
    lineBresenham(front.points[1], front.points[2], color);
    lineBresenham(front.points[2], front.points[3], color);
    lineBresenham(front.points[3], front.points[0], color);
}

void WinInstance::DrawCube(Cube& cube, const sf::Color& color){
    for (auto& face : cube.arr) {
        if (face.n.z() < 0)
            DrawBounds(face, color);
    }
}

void WinInstance::DrawBounds(Cube& cube,const sf::Color& color){
    for (auto& face : cube.arr)
        DrawBounds(face, color);
}

void WinInstance::DrawOnePointProjection(Cube& cube, double r, const sf::Color& color) {
    Point3D<double> center_projection = cube.OnePointTransform(cube.center, r);
    for (auto &face: cube.arr) {
        Point3D<double> face_center = cube.OnePointTransform(face.center, r);
        std::vector<Point3D<double>> points(4);
        for (size_t i = 0; i < face.points.size(); i++)
            points[i] = cube.OnePointTransform(face.points[i], r);

        Point3D<double> n = cross(points[1] - points[0], points[2] - points[1]);
        if (n * (center_projection - face_center) < 0)
            n = -n;
        if (n.z() < 0)
            continue;

        lineBresenham(points[0], points[1], color);
        lineBresenham(points[1], points[2], color);
        lineBresenham(points[2], points[3], color);
        lineBresenham(points[3], points[0], color);
    }
}

