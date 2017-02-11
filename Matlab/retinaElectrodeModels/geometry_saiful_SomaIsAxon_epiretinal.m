function out = model
%
% geometry_saiful_SomaIsAxon_epiretinal.m
%
% Model exported on Apr 11 2015, 11:56 by COMSOL 4.2.0.150.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

model.modelPath('/home/plejaden/Dokumente/Uni/MasterThesis/trunk/src/geometrySources/Saiful_Dokos');

model.name('saiful_SomaIsAxon_epiretinal.mph');

model.param.set('normalDistA', '0.5', 'distance between center and first electrode');
model.param.set('normalDistB', 'normalDistC*0.5');
model.param.set('normalDistC', 'normalDistA/sin(60*2*pi/360)');
model.param.set('electrodeHeight', '10*10^-6');
model.param.set('electrodeZpos', '390*10^-6');

model.modelNode.create('supersimpleRetina');
model.modelNode('supersimpleRetina').name('saiful dokos Model');
model.modelNode('supersimpleRetina').identifier('continuumModel');

model.geom.create('geom1', 3);
model.geom('geom1').feature.create('blk1', 'Block');
model.geom('geom1').feature.create('cyl1', 'Cylinder');
model.geom('geom1').feature.create('blk2', 'Block');
model.geom('geom1').feature.create('blk3', 'Block');
model.geom('geom1').feature.create('blk4', 'Block');
model.geom('geom1').feature.create('blk5', 'Block');
model.geom('geom1').feature.create('blk6', 'Block');
model.geom('geom1').feature.create('blk7', 'Block');
model.geom('geom1').feature.create('blk8', 'Block');
model.geom('geom1').feature.create('blk9', 'Block');
model.geom('geom1').feature.create('blk10', 'Block');
model.geom('geom1').feature.create('blk11', 'Block');
model.geom('geom1').feature.create('blk12', 'Block');
model.geom('geom1').feature.create('blk13', 'Block');
model.geom('geom1').feature.create('blk14', 'Block');
model.geom('geom1').feature.create('wp1', 'WorkPlane');
model.geom('geom1').feature.create('sph1', 'Sphere');
model.geom('geom1').feature.create('cyl2', 'Cylinder');
model.geom('geom1').feature('blk1').name('RGCLayer');
model.geom('geom1').feature('blk1').set('pos', {'0' '0' '270*10^-6'});
model.geom('geom1').feature('blk1').set('size', {'6*10^-3' '4*10^-3' '20*10^-6'});
model.geom('geom1').feature('cyl1').name('RGC_Axon');
model.geom('geom1').feature('cyl1').set('pos', {'1.5*10^-3' '49*10^-6+0.5*10^-3' '285*10^-6'});
model.geom('geom1').feature('cyl1').set('axis', {'0' '1' '0'});
model.geom('geom1').feature('cyl1').set('r', '1*10^-6');
model.geom('geom1').feature('cyl1').set('h', '3*10^-3');
model.geom('geom1').feature('blk2').name('choroid');
model.geom('geom1').feature('blk2').set('pos', {'0' '0' '0'});
model.geom('geom1').feature('blk2').set('size', {'6*10^-3' '4*10^-3' '100*10^-6'});
model.geom('geom1').feature('blk3').name('vitreousLayer');
model.geom('geom1').feature('blk3').set('pos', {'0' '0' '290*10^-6'});
model.geom('geom1').feature('blk3').set('size', {'6*10^-3' '4*10^-3' '200*10^-6'});
model.geom('geom1').feature('blk4').name('electrode_hex_I');
model.geom('geom1').feature('blk4').set('pos', {'1.5 * 10^-3+(-normalDistB*10^-3)' '2* 10^-3+(normalDistA*10^-3)' 'electrodeZpos'});
model.geom('geom1').feature('blk4').set('size', {'10*10^-5' '10*10^-5' 'electrodeHeight'});
model.geom('geom1').feature('blk5').name('electrode_hex_II');
model.geom('geom1').feature('blk5').set('pos', {'1.5 * 10^-3+(+normalDistB*10^-3)' '2* 10^-3+(normalDistA*10^-3)' 'electrodeZpos'});
model.geom('geom1').feature('blk5').set('size', {'10*10^-5' '10*10^-5' 'electrodeHeight'});
model.geom('geom1').feature('blk6').name('electrode_hex_III');
model.geom('geom1').feature('blk6').set('pos', {'1.5 * 10^-3+(+normalDistC*10^-3)' '2* 10^-3+(0*10^-3)' 'electrodeZpos'});
model.geom('geom1').feature('blk6').set('size', {'10*10^-5' '10*10^-5' 'electrodeHeight'});
model.geom('geom1').feature('blk7').name('electrode_hex_IV');
model.geom('geom1').feature('blk7').set('pos', {'1.5 * 10^-3+(+normalDistB*10^-3)' '2* 10^-3+(-normalDistA*10^-3)' 'electrodeZpos'});
model.geom('geom1').feature('blk7').set('size', {'10*10^-5' '10*10^-5' 'electrodeHeight'});
model.geom('geom1').feature('blk8').name('electrode_hex_V');
model.geom('geom1').feature('blk8').set('pos', {'1.5 * 10^-3+(-normalDistB*10^-3)' '2* 10^-3+(-normalDistA*10^-3)' 'electrodeZpos'});
model.geom('geom1').feature('blk8').set('size', {'10*10^-5' '10*10^-5' 'electrodeHeight'});
model.geom('geom1').feature('blk9').name('electrode_hex_VI');
model.geom('geom1').feature('blk9').set('pos', {'1.5 * 10^-3+(-normalDistC*10^-3)' '2* 10^-3+(0*10^-3)' 'electrodeZpos'});
model.geom('geom1').feature('blk9').set('size', {'10*10^-5' '10*10^-5' 'electrodeHeight'});
model.geom('geom1').feature('blk10').name('electrode_center_I');
model.geom('geom1').feature('blk10').set('pos', {'15*10^-4' '20*10^-4' 'electrodeZpos'});
model.geom('geom1').feature('blk10').set('size', {'10*10^-5' '10*10^-5' 'electrodeHeight'});
model.geom('geom1').feature('blk11').name('InnerPlexiformlayer');
model.geom('geom1').feature('blk11').set('pos', {'0' '0' '235*10^-6'});
model.geom('geom1').feature('blk11').set('size', {'6 * 10^-3' '4 * 10^-3' '35* 10^-6'});
model.geom('geom1').feature('blk12').name('nuclearLayers');
model.geom('geom1').feature('blk12').set('pos', {'0' '0' '175 * 10^-6'});
model.geom('geom1').feature('blk12').set('size', {'6 * 10^-3' '4 * 10^-3' '60 * 10^-6'});
model.geom('geom1').feature('blk13').name('subretinalSpace');
model.geom('geom1').feature('blk13').set('pos', {'0' '0' '120 * 10^-6'});
model.geom('geom1').feature('blk13').set('size', {'6 * 10^-3' '4 * 10^-3' '55 * 10^-6'});
model.geom('geom1').feature('blk14').name('retinalPigmentedEpithelium');
model.geom('geom1').feature('blk14').set('pos', {'0' '0' '100 * 10^-6'});
model.geom('geom1').feature('blk14').set('size', {'6 * 10^-3' '4 * 10^-3' '20 * 10^-6'});
model.geom('geom1').feature('wp1').set('quickz', '280*10^-6');
model.geom('geom1').feature('sph1').name('RGC_Soma');
model.geom('geom1').feature('sph1').set('pos', {'1.5*10^-3' '0.5*10^-3' '255*10^-6'});
model.geom('geom1').feature('sph1').set('r', '8*10^-6');
model.geom('geom1').feature('cyl2').name('RGC_IS');
model.geom('geom1').feature('cyl2').set('pos', {'1.5*10^-3' '9*10^-6 + 0.5*10^-3' '266*10^-6'});
model.geom('geom1').feature('cyl2').set('axis', {'0' '0.5' '0.26'});
model.geom('geom1').feature('cyl2').set('r', '1*10^-6');
model.geom('geom1').feature('cyl2').set('h', '40*10^-6');
model.geom('geom1').run;

model.selection.create('sel1', 'Explicit');
model.selection('sel1').geom('geom1', 2);
model.selection('sel1').set([73 74 76 77 78]);
model.selection.create('sel2', 'Explicit');
model.selection('sel2').set([8 9 10 16 17 18]);
model.selection.create('sel3', 'Explicit');
model.selection('sel3').set([1]);
model.selection.create('sel6', 'Explicit');
model.selection('sel6').set([11 12 13 14]);
model.selection.create('sel7', 'Explicit');
model.selection('sel7').set([7]);
model.selection.create('sel8', 'Explicit');
model.selection('sel8').set([6]);
model.selection.create('sel9', 'Explicit');
model.selection('sel9').set([5]);
model.selection.create('sel10', 'Explicit');
model.selection('sel10').set([4]);
model.selection.create('sel11', 'Explicit');
model.selection('sel11').set([3]);
model.selection.create('sel12', 'Explicit');
model.selection('sel12').set([2]);
model.selection.create('sel13', 'Explicit');
model.selection('sel13').set([1]);
model.selection('sel1').name('centerElectrode');
model.selection('sel2').name('hexagonElectrodes');
model.selection('sel3').name('sclera');
model.selection('sel6').name('opticalNerve');
model.selection('sel7').name('vitreous');
model.selection('sel8').name('RGClayer');
model.selection('sel9').name('InnerPlexiformLayer');
model.selection('sel10').name('NuclearLayer');
model.selection('sel11').name('subretinalSpace');
model.selection('sel12').name('pigmentedEpithelium');
model.selection('sel13').name('choroid');

model.variable.create('var1');
model.variable('var1').model('supersimpleRetina');
model.variable('var1').set('retinaConductivity', '1.75');
model.variable('var1').set('a', '1', 'triangle length a');
model.variable('var1').set('b', '0.57735026919', 'triangle length b');
model.variable('var1').set('c', '1.15470053838', 'triangle length c');
model.variable('var1').set('Edist', '10^-3');

model.view('view1').hideObjects.create('hide1');
model.view('view1').hideEntities.create('hide1');

model.material.create('mat1');
model.material('mat1').selection.named('sel8');
model.material.create('mat2');
model.material('mat2').selection.named('sel6');
model.material.create('mat3');
model.material('mat3').propertyGroup.create('linzRes', 'Linearized resistivity');
model.material('mat3').selection.set([8 9 10 15 16 17 18]);
model.material.create('mat5');
model.material('mat5').selection.named('sel7');
model.material.create('mat7');
model.material('mat7').selection.named('sel9');
model.material.create('mat8');
model.material('mat8').selection.named('sel10');
model.material.create('mat9');
model.material('mat9').selection.named('sel11');
model.material.create('mat10');
model.material('mat10').selection.named('sel12');
model.material.create('mat11');
model.material('mat11').selection.named('sel13');

model.physics.create('ec', 'ConductiveMedia', 'geom1');
model.physics('ec').feature.create('gnd1', 'Ground', 2);
model.physics('ec').feature('gnd1').selection.set([30 31 33 34 35 36 37 39 40 41 42 44 46 47 79 80 82 83 84 85 87 88 89 90 91 92 94 95 96]);
model.physics('ec').feature.create('pot1', 'ElectricPotential', 2);
model.physics('ec').feature('pot1').selection.named('sel1');

model.mesh.create('mesh1', 'geom1');
model.mesh('mesh1').feature.create('ftet1', 'FreeTet');

model.result.table.create('evl3', 'Table');

model.variable('var1').name('Variables 1a');

model.view('view1').set('renderwireframe', true);
model.view('view1').set('showgrid', false);
model.view('view1').set('scenelight', 'off');
model.view('view1').hideObjects('hide1').set({'blk1(1)'});

model.material('mat1').name('RGCMaterial');
model.material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat1').propertyGroup('def').set('electricconductivity', {'0.014' '0' '0' '0' '0.014' '0' '0' '0' '0.014'});
model.material('mat1').propertyGroup('def').set('heatcapacity', '820[J/(kg*K)]');
model.material('mat1').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat1').propertyGroup('def').set('emissivity', '0.7');
model.material('mat1').propertyGroup('def').set('density', '2600[kg/m^3]');
model.material('mat1').propertyGroup('def').set('thermalconductivity', {'3[W/(m*K)]' '0' '0' '0' '3[W/(m*K)]' '0' '0' '0' '3[W/(m*K)]'});
model.material('mat2').name('opticalNerveMaterial');
model.material('mat2').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat2').propertyGroup('def').set('electricconductivity', {'1.66' '0' '0' '0' '1.66' '0' '0' '0' '1.66'});
model.material('mat2').propertyGroup('def').set('heatcapacity', '820[J/(kg*K)]');
model.material('mat2').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat2').propertyGroup('def').set('emissivity', '0.7');
model.material('mat2').propertyGroup('def').set('density', '2600[kg/m^3]');
model.material('mat2').propertyGroup('def').set('thermalconductivity', {'3[W/(m*K)]' '0' '0' '0' '3[W/(m*K)]' '0' '0' '0' '3[W/(m*K)]'});
model.material('mat3').name('Copper');
model.material('mat3').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat3').propertyGroup('def').set('electricconductivity', {'5.998e7[S/m]' '0' '0' '0' '5.998e7[S/m]' '0' '0' '0' '5.998e7[S/m]'});
model.material('mat3').propertyGroup('def').set('heatcapacity', '385[J/(kg*K)]');
model.material('mat3').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat3').propertyGroup('def').set('emissivity', '0.5');
model.material('mat3').propertyGroup('def').set('density', '8700[kg/m^3]');
model.material('mat3').propertyGroup('def').set('thermalconductivity', {'400[W/(m*K)]' '0' '0' '0' '400[W/(m*K)]' '0' '0' '0' '400[W/(m*K)]'});
model.material('mat3').propertyGroup('linzRes').set('rho0', '1.72e-8[ohm*m]');
model.material('mat3').propertyGroup('linzRes').set('alpha', '3.9e-3[1/K]');
model.material('mat3').propertyGroup('linzRes').set('Tref', '273.15[K]');
model.material('mat3').propertyGroup('linzRes').addInput('temperature');
model.material('mat5').name('vitreousLayerMaterial');
model.material('mat5').propertyGroup('def').set('electricconductivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat5').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat7').name('InnerPlexiformlayerMaterial');
model.material('mat7').propertyGroup('def').set('electricconductivity', {'0.0549' '0' '0' '0' '0.0549' '0' '0' '0' '0.0549'});
model.material('mat7').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat8').name('nuclearLayersMaterial');
model.material('mat8').propertyGroup('def').set('electricconductivity', {'0.0167' '0' '0' '0' '0.0167' '0' '0' '0' '0.0167'});
model.material('mat8').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat9').name('subretinalSpaceMaterial');
model.material('mat9').propertyGroup('def').set('electricconductivity', {'0.0787' '0' '0' '0' '0.0787' '0' '0' '0' '0.0787'});
model.material('mat9').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat10').name('pigmentedEpitheliumMaterial');
model.material('mat10').propertyGroup('def').set('electricconductivity', {'8.13*10^-4' '0' '0' '0' '8.13*10^-4' '0' '0' '0' '8.13*10^-4'});
model.material('mat10').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat11').name('ChoroidMaterial');
model.material('mat11').propertyGroup('def').set('electricconductivity', {'0.4348' '0' '0' '0' '0.4348' '0' '0' '0' '0.4348'});
model.material('mat11').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});

model.physics('ec').feature('pot1').set('V0', '1');

model.mesh('mesh1').feature('size').set('hauto', 3);
model.mesh('mesh1').run;

model.result.table('evl3').name('Evaluation 3D');
model.result.table('evl3').comments('Interactive 3D values');

model.study.create('std1');
model.study('std1').feature.create('time', 'Transient');
model.study.create('mgmcases');
model.study('mgmcases').feature.create('mgstat', 'Stationary');

model.sol.create('sol1');
model.sol('sol1').study('std1');
model.sol('sol1').attach('std1');
model.sol('sol1').feature.create('st1', 'StudyStep');
model.sol('sol1').feature.create('v1', 'Variables');
model.sol('sol1').feature.create('s1', 'Stationary');
model.sol('sol1').feature('s1').feature.create('fc1', 'FullyCoupled');
model.sol('sol1').feature('s1').feature.create('i1', 'Iterative');
model.sol('sol1').feature('s1').feature('i1').feature.create('mg1', 'Multigrid');
model.sol('sol1').feature('s1').feature.remove('fcDef');

model.result.create('pg1', 'PlotGroup3D');
model.result('pg1').feature.create('slc1', 'Slice');
model.result('pg1').feature.create('slc2', 'Slice');
model.result('pg1').feature.create('slc3', 'Slice');

model.sol('sol1').feature('s1').feature('i1').set('linsolver', 'cg');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('prefun', 'amg');
model.sol('sol1').runAll;

model.result('pg1').name('Electric Potential (ec)');
model.result('pg1').set('title', 'Slice: Electric potential (V) Slice: Electric potential (V) Slice: Electric potential (V) ');
model.result('pg1').set('titleactive', false);
model.result('pg1').feature('slc1').set('quickplane', 'xy');
model.result('pg1').feature('slc1').set('quickznumber', '1');
model.result('pg1').feature('slc2').set('quickxnumber', '1');
model.result('pg1').feature('slc2').set('inheritplot', 'slc1');
model.result('pg1').feature('slc3').set('quickplane', 'zx');
model.result('pg1').feature('slc3').set('quickynumber', '1');
model.result('pg1').feature('slc3').set('inheritplot', 'slc1');

out = model;
