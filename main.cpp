#include <SFML/Graphics.hpp>
#include <iostream>

#include "func.h"

int main() {
  sf::RenderWindow window(sf::VideoMode(800, 600), "Window");

  WinInstance w(window);
  std::vector<Point<int>> points;
  bool Animation = false;
  Point<double> p1(0.0, 0.0);
  Point3D<double> p_3D(150.0,150.0,50.0);
  Cube cube(p_3D,100, 200, 300);
  Cube cube_is(Point3D<double>(350.0,250.0,100.0),100, 100, 100);
  enum Test{
      simple_polygon = 1,
      nonconvex_simple_polygon,
      nonsimple_polygon,
      line,
      hatched_line,
      NZW,
      EO,
      bezier_curve_with_loop,
      bezier_5,
      cyrus_beck,
      catmullromlines,
      zatravka,
      drawbounds,
      drawonepointprojection,
      drawcube,
      isometric,
      animation
  };
  int TEST {16};

  switch(TEST) {
      case simple_polygon:
          points = {Point(150, 100), Point(600, 200), Point(780, 500),
                                         Point(120, 300), Point(150, 100)};
          w.polygon(points, sf::Color::Black);
          if (w.checkConvex(points))
             std::cout << "convex" << std::endl;
          else
             std::cout << "nonconvex" << std::endl;
          if (w.checkPolygonSimplicity(points, p1))
              std::cout << "simple" << std::endl;
          else {
              std::cout << "nonsimple" << std::endl;
              std::cout << p1.x() << " " << p1.y() << std::endl;
          }
          break;
      case nonconvex_simple_polygon:
          points = {Point(100, 100), Point(400, 400), Point(600, 200), Point(780, 500),
                                         Point(300, 570), Point(123, 300), Point(100, 100)};
          w.polygon(points, sf::Color::Black);
          if (w.checkConvex(points))
              std::cout << "convex" << std::endl;
          else
              std::cout << "nonconvex" << std::endl;

          if (w.checkPolygonSimplicity(points, p1))
              std::cout << "simple" << std::endl;
          else {
             std::cout << "nonsimple" << std::endl;
             std::cout << p1.x() << " " << p1.y() << std::endl;
            }
          break;
      case nonsimple_polygon:
          points = {Point(400, 25), Point(700, 550), Point(700, 100),
                                         Point(50, 550), Point(400, 25)};
          w.polygon(points, sf::Color::Black);
          if (w.checkConvex(points))
              std::cout << "convex" << std::endl;
          else
              std::cout << "nonconvex" << std::endl;

          if (w.checkPolygonSimplicity(points, p1))
              std::cout << "simple" << std::endl;
          else {
              std::cout << "nonsimple" << std::endl;
              std::cout << p1.x() << " " << p1.y() << std::endl;
          }
          break;
      case line:
          w.lineBresenham(Point(345, 123), Point(742, 1), sf::Color::Green);
          break;
      case hatched_line:
          w.HatchedLine(Point(70, 30), Point(510, 510), sf::Color::Green);
          break;
      case NZW:
          points = {Point(50, 400), Point(150, 300), Point(250, 400),
                    Point(80, 330), Point(420, 330), Point(173, 500), Point(50, 400)};
          w.fillPolygon(points, methods::NZW, sf::Color::Green);
          break;
      case EO:
          points = {Point(50, 400), Point(150, 300), Point(250, 400),
                    Point(80, 330), Point(420, 330), Point(173, 500), Point(50, 400)};
          w.fillPolygon(points, methods::EO, sf::Color::Green);
          break;

      case bezier_curve_with_loop:
          points = {Point(200, 500), Point(780, 110), Point(20, 100), Point(650, 500)};
          w.curveBezier3(points, sf::Color::Black);
          break;

      case bezier_5:
          points = {Point(100, 100), Point(180, 50), Point(10, 150),
                    Point(280, 100), Point(320, 170),
                    Point(370, 200)};
          w.curveBezier5(points, sf::Color::Black);
          break;


      case cyrus_beck:
          points = {Point(100, 400), Point(600, 400), Point(500, 100), Point(200, 100), Point(100, 400)};
          w.polygon(points, sf::Color::Green);
          {
              Point<int> p1(50, 500), p2(600, 300); // usual
              // Point<int> p1(200, 200), p2(400, 300); // inside
              // Point<int> p1(600, 300), p2(700, 500);  // outside
              w.lineBresenham(p1, p2, sf::Color::Black);
              Point<int> p1_new, p2_new;
              if (w.clipLineCyrusBeck(points, p1, p2, p1_new, p2_new)) {
                  w.lineBresenham(p1_new, p2_new, sf::Color::Red);
              } else {
                  std::cout << "segment is outside the polygon" << std::endl;
              }
          }
          break;

      case catmullromlines: {
          points = { Point(212, 303),Point(375, 374),Point(340, 521), Point(153, 148), Point(343, 41), Point(460,457), Point(337,260) };
          w.CatmullRomLines(points, sf::Color::Green);
          break;
      }

      case zatravka: {
          points = {Point(50, 400), Point(150, 300), Point(250, 400),
                    Point(80, 330), Point(420, 330), Point(173, 500), Point(50, 400)};
          w.fillPolygon(points, methods::NZW, sf::Color::Green);
          w.Zatravka(sf::Color::Green, sf::Color::Red);
          break;
      }

      case drawbounds: {
          cube.rotate(10, cube.GetCenter());
          w.DrawBounds(cube,sf::Color::Red);
          break;
      }

      case drawonepointprojection: {
          cube.rotate(10, cube.GetCenter());
          w.DrawOnePointProjection(cube, 0.001, sf::Color::Red);
          break;
      }

      case drawcube: {
          cube.rotate(10, cube.GetCenter());
          w.DrawCube(cube, sf::Color::Red);
          break;
      }

      case isometric: {
          //w.DrawBounds(cube_is, sf::Color::Green);
          //cube_is.rotate(50, Point3D<double>(0,0,1));
          //cube_is.rotate(45, Point3D<double>(0,1,0));
          //cube_is.rotate(35.26, Point3D<double>(1,0,0));
          //cube_is.rotate(35.26, Point3D<double>(1,0,0));
          cube_is.DoIsometric();
          cube_is.GetCenter().print();
          //cube_is.rotate(100, cube.GetCenter());
          w.DrawBounds(cube_is, sf::Color::Red);
          //w.DrawCube(cube_is, sf::Color::Red);
          break;
      }

      case animation: {
          Animation = true;
      }
  }

enum projection{
      cube_ = 1,
      common,
      onepoint
  };

int PROJECTION {2};
Point3D<double> rotate_vector(1,1,1);


if (Animation){
    switch(PROJECTION){

        case cube_: {
            w.DrawCube(cube, sf::Color::Red);
            w.drawImage();
            while (w.isOpen()) {
                sf::Event event;
                while (w.pollEvent(event)) {
                    if (event.type == sf::Event::KeyPressed) {
                        w.clear(sf::Color::White);
                        w.display();
                        cube.rotate(0.5, rotate_vector);
                        w.DrawCube(cube, sf::Color::Red);
                        w.drawImage();
                    } else if (event.type == sf::Event::Closed)
                        w.close();
                }
            }
            break;
        }

        case common: {
            w.DrawBounds(cube, sf::Color::Red);
            w.drawImage();
            while (w.isOpen()) {
                sf::Event event;
                while (w.pollEvent(event)) {
                    if (event.type == sf::Event::KeyPressed) {
                        w.clear(sf::Color::White);
                        w.display();
                        cube.rotate(0.5, rotate_vector);
                        w.DrawBounds(cube, sf::Color::Red);
                        w.drawImage();
                    } else if (event.type == sf::Event::Closed)
                        w.close();
                }
            }
            break;
        }

        case onepoint: {
            w.DrawOnePointProjection(cube, 0.001, sf::Color::Red);
            w.drawImage();
            while (w.isOpen()) {
                sf::Event event;
                while (w.pollEvent(event)) {
                    if (event.type == sf::Event::KeyPressed) {
                        w.clear(sf::Color::White);
                        w.display();
                        cube.rotate(0.5, rotate_vector);
                        w.DrawOnePointProjection(cube,  0.001, sf::Color::Red);
                        w.drawImage();
                    } else if (event.type == sf::Event::Closed)
                        w.close();
                }
            }
        }
    }

}
else {
    w.saveImage("image.png");
    w.drawImage();

    while (w.isOpen()) {
        sf::Event event;
        while (w.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                w.close();
        }
    }
}
  return 0;
}