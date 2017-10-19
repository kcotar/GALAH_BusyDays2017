FUNCTION get_sme_data, directory=directory, field=field, sobject_id=sobject_id, setup=setup, mode=mode

if not keyword_set(directory) then directory = 'OUTPUT/'
if not keyword_set(field) then field='bmstar'
if not keyword_set(sobject_id) then sobject_id=140708006401203
if not keyword_set(setup) then setup='DR2'
if not keyword_set(mode) then mode='Sp'

restore,directory+field+'_'+strtrim(string(sobject_id), 1)+'_'+setup+'_'+mode+'_SME.out'
sme_out = sme

sme = { $
    field      : field, $
    sobject_id : sobject_id, $
    setup      : setup, $
    mode       : mode, $
    wave       : sme_out.wave, $
    sob        : sme_out.sob, $
    uob        : sme_out.uob, $
    smod       : sme_out.smod $
    }

return, sme

END
