package com.hartwig.hmftools.cobalt.e2e;

import java.nio.file.Path;

import com.hartwig.hmftools.common.genome.tiny.ChromosomeSnippet;
import com.hartwig.hmftools.common.genome.tiny.TinyGenome;

public class CobaltTestGenome extends TinyGenome
{
    public CobaltTestGenome()
    {
        add(new ChromosomeSnippet("chr1", 10_000_000, 1_000_000, snippetPath( "chr1_part_10.txt")));
        add(new ChromosomeSnippet("chr2", 10_000_000, 1_000_000, snippetPath( "chr2_part_10.txt")));
    }

    private Path snippetPath(String snippet)
    {
        return Path.of("src", "test", "resources", "tinygenome", snippet);
    }
}

