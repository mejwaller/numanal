#ifndef CHARTWIDGET1_HPP
#define CHARTWIDGET1_HPP
#include <qwidget.h>
#include <vector>
#include <string>

class QPixmap;

using namespace std;

class Plottable
{    
public:
    Plottable(const vector<float> &x,const vector<float> &y, const string &name, bool line = false)
        :useLine(line),xdata(x),ydata(y),name(name){}
    typedef enum marker
    {
        point = 0,
        cross,
        plus,
        minus,
        asterisk,
        triangle
    } markerType;

    bool lineSeries() const {return useLine;}
    const vector<float>& xData() const {return xdata;}
    const vector<float>& yData() const {return ydata;}
protected:
    bool useLine;
    vector<float> xdata,ydata;
    string name;
};

class ChartWidget1:public QWidget
{
public:
    ChartWidget1(QWidget *parent = 0, const char *name = 0, WFlags f = 0);
    virtual ~ChartWidget1();
    virtual void draw();
    void setNumXTicks(int x) {numxticks = x;}
    void setNumYTicks(int y) {numyticks = y;}
    void setMaxX(float xmax) {maxx = xmax;}
    void setMinX(float xmin) {minx = xmin;}
    void setMaxY(float ymax) {maxy = ymax;}
    void setMinY(float ymin) {miny = ymin;}
    void addPoints(const vector<float> &x, const vector<float> &y,const string &name);
    void addLine(const vector<float> &x, const vector<float> &y,const string &name);
    void drawSeries(QPainter *);
protected:
    void paintEvent( QPaintEvent *);
    QPixmap *buffer;
private:
    void draw(QPainter *);
    float maxx,minx,maxy,miny;
    unsigned int numxticks,numyticks;
    vector<Plottable> series;
    float scale;
};
#endif//CHARTWIDGET1_HPP