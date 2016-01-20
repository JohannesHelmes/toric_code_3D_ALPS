#include <alps/alea/observableset.h>

class measurer {
    measurer(&msmt) {
    }
public:
    void measure();
private:
    ObservableSet& measurements;

}

class thermo_int : measurer {
}

class switching : measurer {
}
