#Spatial SBML Settings
The SpatialUI application uses a custom SBML annotation of the SBML `Model` element, to store the initial setup of a loaded spatial model. The format looks like this: 

      <spatialInfo xmlns="http://fbergmann.github.io/spatial-sbml/settings">
        <update step="0.00001" freq="100"/>
        <items>
          <item sbmlId="p65_P_cytoplasm" palette="black-blue.txt" max="0.001"/>
		  ...
        </items>
      </spatialInfo>

and is described in more detail below.

## The Update Element
The update element holds two attributes: 

* `step` of type `double` describes the initial step size used for the model, and is to be chosen such that the integrator should not fall over immidiately. 
* `freq` of type `int` describes how often the application should attempt to render the scene. 

## The item element
The item element also consists of three elements:

* `sbmlId` of type `string` is the id of a spatial species to be visualized. 
* `palette` of type `string` is the name of a palette to be used. It should be the value of one of the palettes distributed (and present) with `SpatialUI`
* `max` of type `double` is then the maximum value the species will take during the simulation. The values of the palette chosen will be scaled from [0...`max`] and assign `0` the min value of the palette, and `max` accordingly the highest value. 

