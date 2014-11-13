require = function() {};
MykrobeTarget = {};
kTargetSpeciesTB = 1;
MykrobeTarget.species = kTargetSpeciesTB;

describe('Model', function() {
    var model;

	it('can load a model', function(done) {
        model = new Model();
        $(model).on("Model:didLoad", function(event_) {
        	expect(model).not.toEqual(undefined);
            done();
        });
        $(model).on("Model:error", function(event_, error_) {
            console.log('Error: ' + error_.description);
        });
        model.loadFileWithPath('spec/fixtures/tb-2.json');
    });

    it('can determine species and lineage', function() {
        expect(model.speciesModel).toEqual('M. tuberculosis (lineage: European/American)');
    });

    it('can detect mdr', function() {
        expect(model.drugsResistanceModel.mdr).toEqual(true);
    });

    it('can detect xdr', function() {
        expect(model.drugsResistanceModel.xdr).toEqual(true);
    });
});