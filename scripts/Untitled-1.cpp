#include <iostream>
#include <stdexcept>
#include <vector>

struct Never {};
Never never;

union DurationType {
  int i;
  Never n;
};


class Intervention {
  int _duration;
  int _Ngroups;
  int _interventiontimes[]; 

  public:
    //Intervention(int duration, int ngroups, int interventiontimes[]) {
    template <typename T>
    Intervention(int duration, int ngroups, std::vector<T> interventiontimes) {
      if (duration < 1) {
        throw std::invalid_argument("must have a duration of at least 1");
      }
      if (ngroups < 2) {
        throw std::invalid_argument("must have at least 2 groups");
      }
      _duration = duration;
      _Ngroups = ngroups;
      int _interventiontimes[ngroups]; 
      for (int j = 0; j < ngroups; j++) {
        if (interventiontimes[j] <=1) {
          throw std::invalid_argument("all intervention times must be >= 1 (all groups must be initially untreated)");
        } else if (interventiontimes[j] > _duration) {
          _interventiontimes[j] = _duration + 1;
        } else {
          _interventiontimes[j] = interventiontimes[j];
        }
      }
    }

  public:
    int NumberOfGroups() {
      return _Ngroups;
    }

  public:
    int Duration() {
      return _duration;
    }

  public: 
    int Intervened(int group, int t) {
      if (group < 0) {
        throw std::invalid_argument("group cannot be < 0");
      }
      if (group >= _Ngroups) {
        throw std::invalid_argument("group cannot be >= than the number of groups");
      }
      if (t < 0) {
        throw std::invalid_argument("time cannot be < 0");
      }
      if (t > _duration) {
        throw std::invalid_argument("time cannot be > than duration");
      }
      if (t < _interventiontimes[group]) {
        return 0;
      } else {
        return 1;
      }
    }
};

int main() {
  //int newievents[2] = {2, 6};
  std::vector<int> newievents = {2, 6};
  Intervention newi(4, 2, newievents);
  int a = newi.Intervened(0, 3);
  std::cout << "a = "; 
  std::cout << a;
  std::vector<DurationType> asdf = {1, 2, never};
  std::cout << "; asdf = ";
  std::cout << asdf[1];
  return 0;
}


