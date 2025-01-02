package com.hartwig.hmftools.pave.transval;

import java.io.File;
import java.net.URISyntaxException;
import java.net.URL;

import com.hartwig.hmftools.common.genome.refgenome.RefGenomeInterface;

public class TransvalTestBase
{
    public final File ensemblDataDir;
    public final RefGenomeInterface genome = new Chr7Genome();

    public TransvalTestBase()
    {
        URL resourceUrl = getClass().getClassLoader().getResource("ensembl_mini");
        try
        {
            assert resourceUrl != null;
            ensemblDataDir = new File(resourceUrl.toURI());
        }
        catch(URISyntaxException e)
        {
            throw new RuntimeException(e);
        }
    }

    public SingleAminoAcidVariant variant(String definition)
    {
        return new VariationParser(ensemblDataDir).parse(definition);
    }
}
