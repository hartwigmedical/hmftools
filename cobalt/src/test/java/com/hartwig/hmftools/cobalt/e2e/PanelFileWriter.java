package com.hartwig.hmftools.cobalt.e2e;

public class PanelFileWriter extends SectionalData<PanelFileSection>
{
    @Override
    public String header()
    {
        return "chromosome\tposition\trelativeEnrichment";
    }
}
