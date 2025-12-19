/*
  -------------------------------------------------------------------
  
  Copyright (C) 2012-2025, Mahmudul Hasan Anik, Satyajit Roy, and 
  Andrew W. Steiner
  
  This file is part of mcmc-ml.

  Mcmc-ml is free software; you can redistribute it and/or modify it 
  under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  Mcmc-ml is distributed in the hope that it will be useful, but 
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with mcmc-ml. If not, see <http://www.gnu.org/licenses/>.

  -------------------------------------------------------------------
*/

#ifndef MODELS_H
#define MODELS_H

#include "data.h"
#include "settings.h"

namespace mc2ml {

  class model {

  public:

    std::shared_ptr<data> dat;
    std::shared_ptr<settings> set;

    model() {}

    virtual ~model() {}

  };


  class dtr : public model {

  public:

    dtr() {}

    virtual ~dtr() {}

  };


  class dnn : public model {

  public:

    dnn() {}

    virtual ~dnn() {}

  };


  class gp : public model {

  public:

    gp() {}

    virtual ~gp() {}

  };


  class mlpr : public model {

  public:

    mlpr() {}

    virtual ~mlpr() {}

  };

  
  class kde : public model {

  public:

    kde() {}

    virtual ~kde() {}

  };

}

#endif // MODELS_H