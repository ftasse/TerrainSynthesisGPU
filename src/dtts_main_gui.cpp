#include <QApplication>
#include <QtDeclarative/QDeclarativeItem>
#include <QtDeclarative/QDeclarativeView>
#include <QtDeclarative/QDeclarativeComponent>

int main (int argc, char** argv)
{
     QApplication app(argc, argv);
     
     QDeclarativeView view;
     view.setSource(QUrl("qrc:///gui/terrain_synthesis_ui.qml"));
     view.show();
     QObject *object = view.rootObject();
     
     return app.exec();
}