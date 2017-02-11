function out = model
%
% geometry_saiful_SomaIsAxon_epiretinal_300um_a.m
%
% Model exported on Aug 21 2015, 16:03 by COMSOL 4.4.0.150.

import com.comsol.model.*
import com.comsol.model.util.*

model = ModelUtil.create('Model');

%model.modelPath('/home/plejaden/Documents/Uni/MasterThesis/trunk/src/geometrySources/distanceVariationModels_Axon');

%model.name('geometry_saiful_SomaIsAxon_epiretinal_500um_a.mph');

model.comments('Computing the Resistance, voltage and current density of a stimulated retina.');

model.param.set('normalDistA', '0.300/2', 'distance between center and first electrode');
model.param.set('normalDistB', 'normalDistC*0.5');
model.param.set('normalDistC', 'normalDistA/sin(60*2*pi/360)');
model.param.set('electrodeHeight', '10*10^-6');
model.param.set('electrodeEdge', '20*10^-6');
model.param.set('electrodeZpos', '315*10^-6');
model.param.set('neuronYpos', '1.2*10^-3');

model.modelNode.create('comp1');

model.geom.create('geom1', 3);
model.geom('geom1').geomRep('comsol');
model.geom('geom1').selection.create('csel1', 'CumulativeSelection');
model.geom('geom1').selection('csel1').name('electrodes');
model.geom('geom1').selection.create('csel2', 'CumulativeSelection');
model.geom('geom1').selection('csel2').name('electrode_center');
model.geom('geom1').selection.create('csel3', 'CumulativeSelection');
model.geom('geom1').selection('csel3').name('electrode_guards');
model.geom('geom1').selection.create('csel4', 'CumulativeSelection');
model.geom('geom1').selection('csel4').name('Neuron');
model.geom('geom1').feature.create('blk7', 'Block');
model.geom('geom1').feature.create('blk1', 'Block');
model.geom('geom1').feature.create('blk2', 'Block');
model.geom('geom1').feature.create('blk3', 'Block');
model.geom('geom1').feature.create('blk4', 'Block');
model.geom('geom1').feature.create('blk5', 'Block');
model.geom('geom1').feature.create('blk6', 'Block');
model.geom('geom1').feature.create('cyl2', 'Cylinder');
model.geom('geom1').feature.create('cyl3', 'Cylinder');
model.geom('geom1').feature.create('cyl4', 'Cylinder');
model.geom('geom1').feature.create('cyl5', 'Cylinder');
model.geom('geom1').feature.create('cyl6', 'Cylinder');
model.geom('geom1').feature.create('cyl7', 'Cylinder');
model.geom('geom1').feature.create('cyl8', 'Cylinder');
model.geom('geom1').feature.create('wp1', 'WorkPlane');
model.geom('geom1').feature.create('cyl9', 'Cylinder');
model.geom('geom1').feature.create('cyl11', 'Cylinder');
model.geom('geom1').feature.create('cyl10', 'Cylinder');
model.geom('geom1').feature.create('sph1', 'Sphere');
model.geom('geom1').feature.create('igf1', 'IgnoreFaces');
model.geom('geom1').feature('blk7').name('vitreous_layer');
model.geom('geom1').feature('blk7').set('size', {'3*10^-3' '4*10^-3' '25*10^-6'});
model.geom('geom1').feature('blk7').set('pos', {'0' '0' '290*10^-6'});
model.geom('geom1').feature('blk1').name('RGC_layer');
model.geom('geom1').feature('blk1').set('size', {'3*10^-3' '4*10^-3' '20 *10^-6'});
model.geom('geom1').feature('blk1').set('pos', {'0' '0' '270*10^-6'});
model.geom('geom1').feature('blk2').name('Inner_Plexiform_layer');
model.geom('geom1').feature('blk2').set('size', {'3*10^-3' '4*10^-3' '35*10^-6'});
model.geom('geom1').feature('blk2').set('pos', {'0' '0' '235*10^-6'});
model.geom('geom1').feature('blk3').name('Nuclear_layer');
model.geom('geom1').feature('blk3').set('size', {'3*10^-3' '4*10^-3' '60*10^-6'});
model.geom('geom1').feature('blk3').set('pos', {'0' '0' '175*10^-6'});
model.geom('geom1').feature('blk4').name('SubretSpace_layer');
model.geom('geom1').feature('blk4').set('size', {'3*10^-3' '4*10^-3' '55*10^-6'});
model.geom('geom1').feature('blk4').set('pos', {'0' '0' '120*10^-6'});
model.geom('geom1').feature('blk5').name('PigmentEpithelium_layer');
model.geom('geom1').feature('blk5').set('size', {'3*10^-3' '4*10^-3' '20*10^-6'});
model.geom('geom1').feature('blk5').set('pos', {'0' '0' '100*10^-6'});
model.geom('geom1').feature('blk6').name('Choroid_layer');
model.geom('geom1').feature('blk6').set('size', {'3*10^-3' '4*10^-3' '100*10^-6'});
model.geom('geom1').feature('blk6').set('pos', {'0' '0' '0'});
model.geom('geom1').feature('cyl2').name('hexElectrode_center');
model.geom('geom1').feature('cyl2').set('createselection', 'on');
model.geom('geom1').feature('cyl2').set('r', 'electrodeEdge');
model.geom('geom1').feature('cyl2').set('contributeto', 'csel2');
model.geom('geom1').feature('cyl2').set('pos', {'15*10^-4' '20*10^-4' 'electrodeZpos'});
model.geom('geom1').feature('cyl2').set('h', 'electrodeHeight');
model.geom('geom1').feature('cyl3').name('hexElectrode_guard_I');
model.geom('geom1').feature('cyl3').set('createselection', 'on');
model.geom('geom1').feature('cyl3').set('r', 'electrodeEdge');
model.geom('geom1').feature('cyl3').set('contributeto', 'csel3');
model.geom('geom1').feature('cyl3').set('pos', {'1.5 * 10^-3+(-normalDistB*10^-3)' '2* 10^-3+(normalDistA*10^-3)' 'electrodeZpos'});
model.geom('geom1').feature('cyl3').set('h', 'electrodeHeight');
model.geom('geom1').feature('cyl4').name('hexElectrode_guard_II');
model.geom('geom1').feature('cyl4').set('createselection', 'on');
model.geom('geom1').feature('cyl4').set('r', 'electrodeEdge');
model.geom('geom1').feature('cyl4').set('contributeto', 'csel3');
model.geom('geom1').feature('cyl4').set('pos', {'1.5 * 10^-3+(+normalDistB*10^-3)' '2* 10^-3+(normalDistA*10^-3)' 'electrodeZpos'});
model.geom('geom1').feature('cyl4').set('h', 'electrodeHeight');
model.geom('geom1').feature('cyl5').name('hexElectrode_guard_III');
model.geom('geom1').feature('cyl5').set('createselection', 'on');
model.geom('geom1').feature('cyl5').set('r', 'electrodeEdge');
model.geom('geom1').feature('cyl5').set('contributeto', 'csel3');
model.geom('geom1').feature('cyl5').set('pos', {'1.5 * 10^-3+(+normalDistC*10^-3)' '2* 10^-3+(0*10^-3)' 'electrodeZpos'});
model.geom('geom1').feature('cyl5').set('h', 'electrodeHeight');
model.geom('geom1').feature('cyl6').name('hexElectrode_guard_IV');
model.geom('geom1').feature('cyl6').set('createselection', 'on');
model.geom('geom1').feature('cyl6').set('r', 'electrodeEdge');
model.geom('geom1').feature('cyl6').set('contributeto', 'csel3');
model.geom('geom1').feature('cyl6').set('pos', {'1.5 * 10^-3+(+normalDistB*10^-3)' '2* 10^-3+(-normalDistA*10^-3)' 'electrodeZpos'});
model.geom('geom1').feature('cyl6').set('h', 'electrodeHeight');
model.geom('geom1').feature('cyl7').name('hexElectrode_guard_V');
model.geom('geom1').feature('cyl7').set('createselection', 'on');
model.geom('geom1').feature('cyl7').set('r', 'electrodeEdge');
model.geom('geom1').feature('cyl7').set('contributeto', 'csel3');
model.geom('geom1').feature('cyl7').set('pos', {'1.5 * 10^-3+(-normalDistB*10^-3)' '2* 10^-3+(-normalDistA*10^-3)' 'electrodeZpos'});
model.geom('geom1').feature('cyl7').set('h', 'electrodeHeight');
model.geom('geom1').feature('cyl8').name('hexElectrode_guard_VI');
model.geom('geom1').feature('cyl8').set('createselection', 'on');
model.geom('geom1').feature('cyl8').set('r', 'electrodeEdge');
model.geom('geom1').feature('cyl8').set('contributeto', 'csel3');
model.geom('geom1').feature('cyl8').set('pos', {'1.5 * 10^-3+(-normalDistC*10^-3)' '2* 10^-3+(0*10^-3)' 'electrodeZpos'});
model.geom('geom1').feature('cyl8').set('h', 'electrodeHeight');
model.geom('geom1').feature('wp1').set('quickz', 'electrodeZpos - 30*10^-6');
model.geom('geom1').feature('wp1').set('unite', 'on');
model.geom('geom1').feature('cyl9').name('neuron_Axon');
model.geom('geom1').feature('cyl9').set('createselection', 'on');
model.geom('geom1').feature('cyl9').set('r', '1*10^-6');
model.geom('geom1').feature('cyl9').set('contributeto', 'csel4');
model.geom('geom1').feature('cyl9').set('axis', {'0' '1' '0'});
model.geom('geom1').feature('cyl9').set('pos', {'1.5*10^-3' 'neuronYpos + 45*10^-6' '285*10^-6'});
model.geom('geom1').feature('cyl9').set('h', '1.5*10^-3');
model.geom('geom1').feature('cyl11').name('neuron_sodiumBand');
model.geom('geom1').feature('cyl11').set('createselection', 'on');
model.geom('geom1').feature('cyl11').set('r', '1*10^-6');
model.geom('geom1').feature('cyl11').set('contributeto', 'csel4');
model.geom('geom1').feature('cyl11').set('axis', {'0' '1' '0'});
model.geom('geom1').feature('cyl11').set('pos', {'1.5*10^-3' 'neuronYpos' '285*10^-6'});
model.geom('geom1').feature('cyl11').set('h', '45*10^-6');
model.geom('geom1').feature('cyl10').name('neuron_IS');
model.geom('geom1').feature('cyl10').set('createselection', 'on');
model.geom('geom1').feature('cyl10').set('r', '1*10^-6');
model.geom('geom1').feature('cyl10').set('contributeto', 'csel4');
model.geom('geom1').feature('cyl10').set('axis', {'0' '0.5' '0.26'});
model.geom('geom1').feature('cyl10').set('pos', {'1.5*10^-3' 'neuronYpos + -35*10^-6' '266*10^-6'});
model.geom('geom1').feature('cyl10').set('h', '40*10^-6');
model.geom('geom1').feature('sph1').name('neuron_soma');
model.geom('geom1').feature('sph1').set('r', '8*10^-6');
model.geom('geom1').feature('sph1').set('pos', {'1.5*10^-3' 'neuronYpos + -40*10^-6' '260*10^-6'});
model.geom('geom1').feature('sph1').set('contributeto', 'csel4');
model.geom('geom1').feature('sph1').set('createselection', 'on');
model.geom('geom1').feature('igf1').selection('input').named('csel4');
model.geom('geom1').run;

model.material.create('mat1');
model.material('mat1').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.material('mat1').propertyGroup.create('linzRes', 'Linearized resistivity');
model.material('mat1').selection.set([6 7 8 9 10 11 12 13 14]);
model.material.create('mat2');
model.material('mat2').propertyGroup('def').func.create('dL_solid_1', 'Piecewise');
model.material('mat2').propertyGroup('def').func.create('k', 'Piecewise');
model.material('mat2').propertyGroup('def').func.create('res_solid_1', 'Piecewise');
model.material('mat2').propertyGroup('def').func.create('epsilon', 'Piecewise');
model.material('mat2').propertyGroup('def').func.create('alpha_solid_1', 'Piecewise');
model.material('mat2').propertyGroup('def').func.create('C_solid_1', 'Piecewise');
model.material('mat2').propertyGroup('def').func.create('HC_solid_1', 'Piecewise');
model.material('mat2').propertyGroup('def').func.create('mu', 'Piecewise');
model.material('mat2').propertyGroup('def').func.create('sigma_solid_1', 'Piecewise');
model.material('mat2').propertyGroup('def').func.create('rho_solid_1', 'Piecewise');
model.material('mat2').propertyGroup('def').func.create('VP_solid_1', 'Piecewise');
model.material('mat2').propertyGroup('def').func.create('kappa', 'Piecewise');
model.material('mat2').propertyGroup.create('Enu', 'Young''s modulus and Poisson''s ratio');
model.material('mat2').propertyGroup('Enu').func.create('E', 'Piecewise');
model.material('mat2').propertyGroup('Enu').func.create('nu', 'Piecewise');
model.material('mat2').selection.set([8 9 10 11 12 13 14]);
model.material.create('mat9');
model.material('mat9').selection.set([7]);
model.material.create('mat3');
model.material('mat3').selection.set([6]);
model.material.create('mat4');
model.material('mat4').selection.set([5]);
model.material.create('mat5');
model.material('mat5').selection.set([4]);
model.material.create('mat6');
model.material('mat6').selection.set([3]);
model.material.create('mat7');
model.material('mat7').selection.set([2]);
model.material.create('mat8');
model.material('mat8').selection.set([1]);

model.physics.create('ec', 'ConductiveMedia', 'geom1');
model.physics('ec').feature.create('gnd1', 'Ground', 2);
model.physics('ec').feature('gnd1').selection.set([33 39 43 57 61 69]);
model.physics('ec').feature.create('term1', 'Terminal', 2);
model.physics('ec').feature('term1').selection.set([51]);

model.mesh.create('mesh1', 'geom1');
model.mesh('mesh1').feature.create('edg1', 'Edge');
model.mesh('mesh1').feature.create('ftet3', 'FreeTet');
model.mesh('mesh1').feature.create('ftet4', 'FreeTet');
model.mesh('mesh1').feature.create('ftet5', 'FreeTet');
model.mesh('mesh1').feature.create('ftet1', 'FreeTet');
model.mesh('mesh1').feature('edg1').selection.named('geom1_csel1_edg');
model.mesh('mesh1').feature('edg1').feature.create('dis1', 'Distribution');
model.mesh('mesh1').feature('edg1').feature.create('dis2', 'Distribution');
model.mesh('mesh1').feature('edg1').feature('dis1').selection.named('geom1_csel3_edg');
model.mesh('mesh1').feature('edg1').feature('dis2').selection.named('geom1_csel2_edg');
model.mesh('mesh1').feature('ftet3').selection.named('geom1_csel2_dom');
model.mesh('mesh1').feature('ftet3').feature.create('size1', 'Size');
model.mesh('mesh1').feature('ftet3').feature('size1').selection.named('geom1_csel2_dom');
model.mesh('mesh1').feature('ftet4').selection.named('geom1_csel3_dom');
model.mesh('mesh1').feature('ftet4').feature.create('size1', 'Size');
model.mesh('mesh1').feature('ftet4').feature('size1').selection.named('geom1_csel3_dom');
model.mesh('mesh1').feature('ftet5').selection.geom('geom1', 3);
model.mesh('mesh1').feature('ftet5').selection.set([6 7]);
model.mesh('mesh1').feature('ftet5').feature.create('size1', 'Size');
model.mesh('mesh1').feature('ftet1').feature.create('size1', 'Size');

model.view('view1').set('renderwireframe', true);
model.view('view1').set('scenelight', 'off');
model.view('view2').axis.set('xmin', '-0.2899113893508911');
model.view('view2').axis.set('ymin', '0');
model.view('view2').axis.set('xmax', '1.2899112701416016');

model.material('mat1').name('Copper');
model.material('mat1').propertyGroup('def').set('relpermeability', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat1').propertyGroup('def').set('electricconductivity', {'5.998e7[S/m]' '0' '0' '0' '5.998e7[S/m]' '0' '0' '0' '5.998e7[S/m]'});
model.material('mat1').propertyGroup('def').set('thermalexpansioncoefficient', {'17e-6[1/K]' '0' '0' '0' '17e-6[1/K]' '0' '0' '0' '17e-6[1/K]'});
model.material('mat1').propertyGroup('def').set('heatcapacity', '385[J/(kg*K)]');
model.material('mat1').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat1').propertyGroup('def').set('density', '8700[kg/m^3]');
model.material('mat1').propertyGroup('def').set('thermalconductivity', {'400[W/(m*K)]' '0' '0' '0' '400[W/(m*K)]' '0' '0' '0' '400[W/(m*K)]'});
model.material('mat1').propertyGroup('Enu').set('youngsmodulus', '110e9[Pa]');
model.material('mat1').propertyGroup('Enu').set('poissonsratio', '0.35');
model.material('mat1').propertyGroup('linzRes').set('rho0', '');
model.material('mat1').propertyGroup('linzRes').set('alpha', '');
model.material('mat1').propertyGroup('linzRes').set('Tref', '');
model.material('mat1').propertyGroup('linzRes').set('rho0', '1.72e-8[ohm*m]');
model.material('mat1').propertyGroup('linzRes').set('alpha', '0.0039[1/K]');
model.material('mat1').propertyGroup('linzRes').set('Tref', '298[K]');
model.material('mat1').propertyGroup('linzRes').addInput('temperature');
model.material('mat2').name('Platinum [solid]');
model.material('mat2').propertyGroup('def').func('dL_solid_1').set('pieces', {'10.0' '70.0' '-0.00192855+1.134373E-6*T^1-9.005848E-8*T^2+3.271584E-9*T^3-3.362442E-11*T^4+1.238208E-13*T^5'; '70.0' '280.0' '-0.001895358-2.225109E-6*T^1+7.277336E-8*T^2-2.321766E-10*T^3+2.913888E-13*T^4'; '280.0' '1973.0' '-0.002617286+8.777313E-6*T^1+3.786334E-10*T^2+5.373731E-13*T^3-6.636722E-17*T^4'});
model.material('mat2').propertyGroup('def').func('dL_solid_1').set('arg', 'T');
model.material('mat2').propertyGroup('def').func('k').set('pieces', {'0.0' '13.0' '209.6991*T^1+18.16709*T^2-4.678988*T^3+0.2278691*T^4-0.002748296*T^5'; '13.0' '50.0' '2978.797-201.1759*T^1+3.362113*T^2+0.0713319*T^3-0.002761393*T^4+2.280531E-5*T^5'; '50.0' '100.0' '1212.843-69.25658*T^1+1.763533*T^2-0.0228956*T^3+1.494028E-4*T^4-3.889706E-7*T^5'; '100.0' '285.0' '123.3886-1.066855*T^1+0.009646914*T^2-4.536814E-5*T^3+1.072687E-7*T^4-1.004522E-10*T^5'; '285.0' '2045.0' '73.99627-0.01557887*T^1+2.646931E-5*T^2-6.133801E-9*T^3'});
model.material('mat2').propertyGroup('def').func('k').set('arg', 'T');
model.material('mat2').propertyGroup('def').func('res_solid_1').set('pieces', {'14.0' '47.0' '-5.497611E-10+1.415797E-10*T^1-1.305787E-11*T^2+5.545444E-13*T^3-7.41428E-15*T^4+3.513902E-17*T^5'; '47.0' '160.0' '1.845544E-9-3.44008E-10*T^1+1.431636E-11*T^2-1.250757E-13*T^3+5.330375E-16*T^4-8.95938E-19*T^5'; '160.0' '600.0' '-1.927892E-8+5.233699E-10*T^1-4.107885E-13*T^2+6.694129E-16*T^3-4.447775E-19*T^4'; '600.0' '2000.0' '-4.843579E-8+5.552497E-10*T^1-1.600249E-13*T^2+2.814022E-17*T^3'});
model.material('mat2').propertyGroup('def').func('res_solid_1').set('arg', 'T');
model.material('mat2').propertyGroup('def').func('epsilon').set('pieces', {'1000.0' '2000.0' '0.1248438+6.688811E-5*T^1+5.827506E-10*T^2'});
model.material('mat2').propertyGroup('def').func('epsilon').set('arg', 'T');
model.material('mat2').propertyGroup('def').func('alpha_solid_1').set('pieces', {'10.0' '70.0' '6.594988E-6+1.839984E-8*T^1+3.921862E-10*T^2-1.081641E-11*T^3+1.007161E-13*T^4-3.398394E-16*T^5'; '70.0' '230.0' '6.867919E-6+1.945355E-8*T^1-6.698432E-11*T^2+1.199734E-13*T^3-1.069967E-16*T^4'; '230.0' '1973.0' '8.801519E-6+4.097477E-10*T^1+1.248065E-12*T^2-7.133932E-16*T^3+1.689741E-19*T^4'});
model.material('mat2').propertyGroup('def').func('alpha_solid_1').set('arg', 'T');
model.material('mat2').propertyGroup('def').func('C_solid_1').set('pieces', {'0.0' '19.0' '0.03281349*T^1+0.001129466*T^2+3.449445E-4*T^3+5.174165E-5*T^4-1.325633E-6*T^5'; '19.0' '119.0' '10.30393-1.986516*T^1+0.1283953*T^2-0.002010741*T^3+1.359791E-5*T^4-3.445457E-8*T^5'; '119.0' '290.0' '0.4467027+1.721765*T^1-0.009418853*T^2+2.453936E-5*T^3-2.455881E-8*T^4'; '290.0' '2000.0' '122.2187+0.03986346*T^1-1.836174E-5*T^2+7.556773E-9*T^3'});
model.material('mat2').propertyGroup('def').func('C_solid_1').set('arg', 'T');
model.material('mat2').propertyGroup('def').func('HC_solid_1').set('pieces', {'0.0' '19.0' '0.006401583*T^1+2.203475E-4*T^2+6.72952E-5*T^3+1.009428E-5*T^4-2.586177E-7*T^5'; '19.0' '119.0' '2.018562-0.3875494*T^1+0.02504865*T^2-3.922754E-4*T^3+2.652817E-6*T^4-6.721742E-9*T^5'; '119.0' '290.0' '0.08714724+0.3358992*T^1-0.001837525*T^2+4.787383E-6*T^3-4.791178E-9*T^4'; '290.0' '2000.0' '23.84364+0.007776964*T^1-3.582192E-6*T^2+1.47425E-9*T^3'});
model.material('mat2').propertyGroup('def').func('HC_solid_1').set('arg', 'T');
model.material('mat2').propertyGroup('def').func('mu').set('pieces', {'293.0' '1480.0' '6.55395E10-1.15E7*T^1+3.035766E-9*T^2'});
model.material('mat2').propertyGroup('def').func('mu').set('arg', 'T');
model.material('mat2').propertyGroup('def').func('sigma_solid_1').set('pieces', {'14.0' '47.0' '1/(3.513902E-17*T^5-7.414280E-15*T^4+5.545444E-13*T^3-1.305787E-11*T^2+1.415797E-10*T-5.497611E-10)'; '47.0' '160.0' '1/(-8.959380E-19*T^5+5.330375E-16*T^4-1.250757E-13*T^3+1.431636E-11*T^2-3.440080E-10*T+1.845544E-09)'; '160.0' '600.0' '1/(-4.447775E-19*T^4+6.694129E-16*T^3-4.107885E-13*T^2+5.233699E-10*T-1.927892E-08)'; '600.0' '2000.0' '1/(2.814022E-17*T^3-1.600249E-13*T^2+5.552497E-10*T-4.843579E-08)'});
model.material('mat2').propertyGroup('def').func('sigma_solid_1').set('arg', 'T');
model.material('mat2').propertyGroup('def').func('rho_solid_1').set('pieces', {'10.0' '70.0' '21512.45-0.07336997*T^1+0.005824972*T^2-2.116082E-4*T^3+2.17523E-6*T^4-8.010422E-9*T^5'; '70.0' '280.0' '21510.32+0.1433013*T^1-0.00470128*T^2+1.502256E-5*T^3-1.883994E-8*T^4'; '280.0' '1973.0' '21557.19-0.5675783*T^1-1.7525E-5*T^2-3.171806E-8*T^3+4.698968E-12*T^4'});
model.material('mat2').propertyGroup('def').func('rho_solid_1').set('arg', 'T');
model.material('mat2').propertyGroup('def').func('VP_solid_1').set('pieces', {'293.0' '2041.0' '(exp((-2.938700e+004/T+1.103900e+000*log10(T)+7.762810e+000-4.527000e-001/T^3)*log(10.0)))*1.333200e+002'});
model.material('mat2').propertyGroup('def').func('VP_solid_1').set('arg', 'T');
model.material('mat2').propertyGroup('def').func('kappa').set('pieces', {'293.0' '1480.0' '1.99731E11-6.619246E7*T^1+4690.398*T^2'});
model.material('mat2').propertyGroup('def').func('kappa').set('arg', 'T');
model.material('mat2').propertyGroup('def').set('dL', '(dL_solid_1(T[1/K])-dL_solid_1(Tempref[1/K]))/(1+dL_solid_1(Tempref[1/K]))');
model.material('mat2').propertyGroup('def').set('thermalconductivity', {'k(T[1/K])[W/(m*K)]' '0' '0' '0' 'k(T[1/K])[W/(m*K)]' '0' '0' '0' 'k(T[1/K])[W/(m*K)]'});
model.material('mat2').propertyGroup('def').set('resistivity', {'res_solid_1(T[1/K])[ohm*m]' '0' '0' '0' 'res_solid_1(T[1/K])[ohm*m]' '0' '0' '0' 'res_solid_1(T[1/K])[ohm*m]'});
model.material('mat2').propertyGroup('def').set('emissivity', 'epsilon(T[1/K])');
model.material('mat2').propertyGroup('def').set('thermalexpansioncoefficient', {'(alpha_solid_1(T[1/K])[1/K]+(Tempref-293[K])*if(abs(T-Tempref)>1e-3,(alpha_solid_1(T[1/K])[1/K]-alpha_solid_1(Tempref[1/K])[1/K])/(T-Tempref),d(alpha_solid_1(T[1/K]),T)[1/K]))/(1+alpha_solid_1(Tempref[1/K])[1/K]*(Tempref-293[K]))' '0' '0' '0' '(alpha_solid_1(T[1/K])[1/K]+(Tempref-293[K])*if(abs(T-Tempref)>1e-3,(alpha_solid_1(T[1/K])[1/K]-alpha_solid_1(Tempref[1/K])[1/K])/(T-Tempref),d(alpha_solid_1(T[1/K]),T)[1/K]))/(1+alpha_solid_1(Tempref[1/K])[1/K]*(Tempref-293[K]))' '0' '0' '0' '(alpha_solid_1(T[1/K])[1/K]+(Tempref-293[K])*if(abs(T-Tempref)>1e-3,(alpha_solid_1(T[1/K])[1/K]-alpha_solid_1(Tempref[1/K])[1/K])/(T-Tempref),d(alpha_solid_1(T[1/K]),T)[1/K]))/(1+alpha_solid_1(Tempref[1/K])[1/K]*(Tempref-293[K]))'});
model.material('mat2').propertyGroup('def').set('heatcapacity', 'C_solid_1(T[1/K])[J/(kg*K)]');
model.material('mat2').propertyGroup('def').set('HC', 'HC_solid_1(T[1/K])[J/(mol*K)]');
model.material('mat2').propertyGroup('def').set('mu', 'mu(T[1/K])[Pa]');
model.material('mat2').propertyGroup('def').set('electricconductivity', {'sigma_solid_1(T[1/K])[S/m]' '0' '0' '0' 'sigma_solid_1(T[1/K])[S/m]' '0' '0' '0' 'sigma_solid_1(T[1/K])[S/m]'});
model.material('mat2').propertyGroup('def').set('density', 'rho_solid_1(T[1/K])[kg/m^3]');
model.material('mat2').propertyGroup('def').set('VP', 'VP_solid_1(T[1/K])[Pa]');
model.material('mat2').propertyGroup('def').set('kappa', 'kappa(T[1/K])[Pa]');
model.material('mat2').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat2').propertyGroup('def').addInput('temperature');
model.material('mat2').propertyGroup('def').addInput('strainreferencetemperature');
model.material('mat2').propertyGroup('Enu').func('E').set('pieces', {'93.0' '293.0' '1.667479E11-4.099535E7*T^1+54115.62*T^2-51.62897*T^3'; '293.0' '1480.0' '1.68E11-3.38E7*T^1'});
model.material('mat2').propertyGroup('Enu').func('E').set('arg', 'T');
model.material('mat2').propertyGroup('Enu').func('nu').set('pieces', {'293.0' '1480.0' '0.3516936-1.897311E-5*T^1-5.685048E-9*T^2'});
model.material('mat2').propertyGroup('Enu').func('nu').set('arg', 'T');
model.material('mat2').propertyGroup('Enu').set('youngsmodulus', 'E(T[1/K])[Pa]');
model.material('mat2').propertyGroup('Enu').set('poissonsratio', 'nu(T[1/K])');
model.material('mat2').propertyGroup('Enu').addInput('temperature');
model.material('mat9').name('Vitreous_material');
model.material('mat9').propertyGroup('def').set('electricconductivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat9').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat3').name('RGC_material');
model.material('mat3').propertyGroup('def').set('electricconductivity', {'1/70.4' '0' '0' '0' '1/70.4' '0' '0' '0' '1/70.4'});
model.material('mat3').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat4').name('Inner_Plexiform_material');
model.material('mat4').propertyGroup('def').set('electricconductivity', {'1/18.2' '0' '0' '0' '1/18.2' '0' '0' '0' '1/18.2'});
model.material('mat4').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat5').name('Nuclear_material');
model.material('mat5').propertyGroup('def').set('electricconductivity', {'1/60' '0' '0' '0' '1/60' '0' '0' '0' '1/60'});
model.material('mat5').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat6').name('Subretinal_material');
model.material('mat6').propertyGroup('def').set('electricconductivity', {'1/12.7' '0' '0' '0' '1/12.7' '0' '0' '0' '1/12.7'});
model.material('mat6').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat7').name('PigmentEpithelium_material');
model.material('mat7').propertyGroup('def').set('electricconductivity', {'1/4370' '0' '0' '0' '1/4370' '0' '0' '0' '1/4370'});
model.material('mat7').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});
model.material('mat8').name('Choroid_material');
model.material('mat8').propertyGroup('def').set('electricconductivity', {'1/2.3' '0' '0' '0' '1/2.3' '0' '0' '0' '1/2.3'});
model.material('mat8').propertyGroup('def').set('relpermittivity', {'1' '0' '0' '0' '1' '0' '0' '0' '1'});

model.physics('ec').feature('term1').set('I0', '1');
model.physics('ec').feature('term1').set('TerminalType', 'Voltage');

model.mesh('mesh1').feature('size').set('hauto', 3);
model.mesh('mesh1').feature('edg1').feature('dis1').set('numelem', '50');
model.mesh('mesh1').feature('edg1').feature('dis2').set('numelem', '50');
model.mesh('mesh1').feature('ftet3').name('Electrdrode_center_tetrahedral');
model.mesh('mesh1').feature('ftet3').feature('size1').set('hauto', 2);
model.mesh('mesh1').feature('ftet4').name('Electrdrode_guard_tetrahedral');
model.mesh('mesh1').feature('ftet4').feature('size1').set('hauto', 2);
model.mesh('mesh1').feature('ftet5').feature('size1').set('hauto', 2);
model.mesh('mesh1').feature('ftet1').feature('size1').set('hauto', 4);
model.mesh('mesh1').run;

model.study.create('std1');
model.study('std1').feature.create('stat', 'Stationary');

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

model.study('std1').feature('stat').set('initstudyhide', 'on');
model.study('std1').feature('stat').set('initsolhide', 'on');
model.study('std1').feature('stat').set('notstudyhide', 'on');
model.study('std1').feature('stat').set('notsolhide', 'on');

model.result.numerical.create('gev1', 'EvalGlobal');
model.result.numerical.create('int1', 'IntSurface');
model.result.numerical.create('int4', 'IntSurface');
model.result.numerical.create('int2', 'IntSurface');
model.result.numerical.create('int3', 'IntSurface');
model.result.numerical('gev1').set('probetag', 'none');
model.result.numerical('int1').selection.set([51]);
model.result.numerical('int1').set('probetag', 'none');
model.result.numerical('int4').selection.set([33 39 43 57 61 69]);
model.result.numerical('int4').set('probetag', 'none');
model.result.numerical('int2').selection.set([43]);
model.result.numerical('int2').set('probetag', 'none');
model.result.numerical('int3').selection.set([16 17 18 28 30 31 34 35 36 37 40 41 44 45 46 47 48 49 52 53 54 55 58 59 62 63 64 65 66 67 70 71 77]);
model.result.numerical('int3').set('probetag', 'none');
model.result.create('pg1', 'PlotGroup3D');
model.result('pg1').feature.create('mslc1', 'Multislice');
model.result('pg1').feature.create('mslc2', 'Multislice');

model.sol('sol1').attach('std1');
model.sol('sol1').feature('st1').name('Compile Equations: Stationary');
model.sol('sol1').feature('st1').set('studystep', 'stat');
model.sol('sol1').feature('st1').set('splitcomplex', true);
model.sol('sol1').feature('v1').set('control', 'stat');
model.sol('sol1').feature('s1').set('keeplog', true);
model.sol('sol1').feature('s1').set('stol', '10^-6');
model.sol('sol1').feature('s1').set('control', 'stat');
model.sol('sol1').feature('s1').feature('fc1').set('maxiter', '100');
model.sol('sol1').feature('s1').feature('i1').set('linsolver', 'cg');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('prefun', 'amg');
model.sol('sol1').feature('s1').feature('i1').feature('mg1').set('iter', '1');
model.sol('sol1').runAll;

model.result.numerical('gev1').set('expr', 'ec.R11');
model.result.numerical('gev1').set('descr', 'Resistance');
model.result.numerical('int1').name('FluxSource');
model.result.numerical('int1').set('method', 'integration');
model.result.numerical('int1').set('unit', 'nA');
model.result.numerical('int1').set('descr', 'Flux_soruce');
model.result.numerical('int1').set('expr', 'nx*ec.Jx+ny*ec.Jy+nz*ec.Jz');
model.result.numerical('int1').set('descractive', true);
model.result.numerical('int1').set('descr', 'Flux_soruce');
model.result.numerical('int4').name('FluxDrain_all');
model.result.numerical('int4').set('unit', 'nA');
model.result.numerical('int4').set('descr', 'Flux_6drain');
model.result.numerical('int4').set('expr', 'nx*ec.Jx+ny*ec.Jy+nz*ec.Jz');
model.result.numerical('int4').set('descractive', true);
model.result.numerical('int4').set('descr', 'Flux_6drain');
model.result.numerical('int2').name('FluxDrain_single');
model.result.numerical('int2').set('unit', 'nA');
model.result.numerical('int2').set('descr', 'Flux_1drain');
model.result.numerical('int2').set('expr', 'nx*ec.Jx+ny*ec.Jy+nz*ec.Jz');
model.result.numerical('int2').set('descractive', true);
model.result.numerical('int2').set('descr', 'Flux_1drain');
model.result.numerical('int3').name('FluxBorder');
model.result.numerical('int3').set('unit', 'nA');
model.result.numerical('int3').set('descr', 'Flux_IsolatedBoarders');
model.result.numerical('int3').set('expr', 'nx*ec.Jx+ny*ec.Jy+nz*ec.Jz');
model.result.numerical('int3').set('descractive', true);
model.result.numerical('int3').set('descr', 'Flux_IsolatedBoarders');
model.result('pg1').name('Electric Potential (ec)');
model.result('pg1').set('frametype', 'spatial');
model.result('pg1').feature('mslc2').active(false);
model.result('pg1').feature('mslc2').set('descr', 'ec.Jx+ec.Jy+ec.Jz');
model.result('pg1').feature('mslc2').set('unit', 'A/m^2');
model.result('pg1').feature('mslc2').set('expr', 'ec.Jx+ec.Jy+ec.Jz');

out = model;
