package com.hartwig.hmftools.iclusion.data;

public class IclusionTrial {

    private final String title;
    private final String acronym;
    private final String eudra;
    private final String nct;
    private final String ipn;
    private final String ccmo;

    public IclusionTrial(final String title, final String acronym, final String eudra, final String nct, final String ipn,
            final String ccmo) {
        this.title = title;
        this.acronym = acronym;
        this.eudra = eudra;
        this.nct = nct;
        this.ipn = ipn;
        this.ccmo = ccmo;
    }
}
