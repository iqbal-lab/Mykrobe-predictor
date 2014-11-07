var targetName = require('./package.json').targetName,
	MykrobeTarget = MykrobeTarget || {},
	kTargetTypePredictor = 0,
	kTargetSpeciesSAureus = 0,
	kTargetSpeciesTB = 1;

MykrobeTarget.targetName = targetName;

if ( targetName === 'predictor-s-aureus' ) {
	MykrobeTarget.type = kTargetTypePredictor;
	MykrobeTarget.species = kTargetSpeciesSAureus;
}
else if ( targetName === 'predictor-tb' ) {
	MykrobeTarget.type = kTargetTypePredictor;
	MykrobeTarget.species = kTargetSpeciesTB;
}
else {
	// default;
	MykrobeTarget.type = kTargetTypePredictor;
	MykrobeTarget.species = kTargetSpeciesSAureus;
}
