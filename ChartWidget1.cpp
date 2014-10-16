#include "ChartWidget1.hpp"

#include <math.h>
#include <qpainter.h>
#include <qpixmap.h>
#include <iostream>

using namespace std;

ChartWidget1::ChartWidget1(QWidget *parent, const char *name, WFlags f)
    :QWidget(parent,name,f),buffer(0),maxx(100.0),minx(0.0),maxy(100.0),miny(0.0),numxticks(1),numyticks(1),scale(100000.0)
{
    setBackgroundMode( QWidget::PaletteBase );       

}

ChartWidget1::~ChartWidget1()
{
    if(buffer)
    {
        delete buffer;
    }
}

void ChartWidget1::paintEvent(QPaintEvent *)
{
    if (buffer && width() == buffer->width() && height() == buffer->height())
    {
        bitBlt(this,0,0,buffer);
    }
    else
    {
        draw();        
    }
}

void ChartWidget1::draw()
{
    if(buffer)
    {
        delete buffer;
    }
    buffer = new QPixmap(width(),height());
    buffer->fill();
    QPainter bufferpaint(buffer);
    draw(&bufferpaint);
    bitBlt(this,0,0,buffer);
}

void ChartWidget1::draw(QPainter *paint)
{    
    QRect v;
    
    //lets leave 10% margin either side...
    float scalefac = 0.8;
    paint->save();   

    //we want to draw so that the graph origin is 10% in eitehr side,
    //so we want coords to go from minx - 10% of width (10% of width = 0.1*maxx-minx) to maxx + 10% of width
    //so window width = maxx-minx*1.2,height = maxy-miny*1.2
    //and similarly for height (mutliplying by scale - but I assume paint->scale can do that?:
    paint->setWindow(scale*(minx-0.1*(maxx-minx)),scale*(miny-0.1*(maxy-miny)),scale*((maxx-minx)*1.2),scale*((maxy-miny)*1.2));
    
    QWMatrix mat1;        
    /*paint->*/mat1.scale(1,-1);
    /*paint->*/mat1.translate(0,-scale*(maxy-miny));
    paint->setWorldMatrix(mat1);    
    //Get the viewport and set it up
    /*v = paint->viewport();
    int d = QMIN( v.width(), v.height() );        

    paint->setViewport( v.left() + (v.width()-d)/2,
        v.top() + (v.height()-d)/2, d, d );    
    v = paint->viewport();   
    */

    paint->drawLine(scale*minx,scale*miny,scale*minx,scale*maxy);//y axis
    paint->drawLine(scale*minx,scale*miny,scale*maxx,scale*miny);//xaxis
    //now draw ticks:
    for(unsigned int i = 0; i< numyticks;i++)
    {
        paint->drawLine(scale*minx,(i+1)*scale*(maxy-miny)/numyticks,-scale*(maxx-minx)*0.01,(i+1)*scale*(maxy-miny)/numyticks);
    }
    for(unsigned int j = 0; j<numxticks;j++)
    {
        paint->drawLine((j+1)*scale*(maxx-minx)/numxticks,scale*miny,(j+1)*scale*(maxx-minx)/numxticks,-scale*(maxy-miny)*0.01);
    }

    drawSeries(paint);
    //paint->drawText(minx,miny,"(minx,miny)");
    //paint->drawText(maxx,maxy,"(maxx,maxy)");
    paint->restore();
    
}


void ChartWidget1::addPoints(const vector<float>& x, const vector<float>& y, const string& name)
{
    Plottable PointPlot(x,y,name);
    series.push_back(PointPlot);

}

void ChartWidget1::addLine(const vector<float>& x, const vector<float>& y, const string& name)
{
    Plottable LinePlot(x,y,name,true);
    series.push_back(LinePlot);
}

void ChartWidget1::drawSeries(QPainter *paint)
{
    vector<float> x,y;
    if(!series.empty())
    {
        paint->setPen(QPen(black,0));
        for(int i = 0; i<series.size(); i++)
        {
            x = series[i].xData();
            y = series[i].yData();
            if(series[i].lineSeries())
            {                
                for(unsigned int j = 0; j < x.size()-1; j++)
                {
                    paint->drawLine(x[j]*scale,y[j]*scale,x[j+1]*scale,y[j+1]*scale);
                }
            }
            else
            {                
                for(unsigned int j = 0; j < x.size(); j++)
                {
                    //paint->drawPoint(x[j]*scale,y[j]*scale);
                    //paint->drawEllipse(x[j]*scale,y[j]*scale,scale/100,scale/100);                    
                    //paint->setPen(QPen(red,50));
                    paint->drawEllipse(x[j]*scale-(scale/100.0)/2.0,y[j]*scale-(scale/100.0)/2.0,scale/100,scale/100);
                }
            }
        }
    }
}







