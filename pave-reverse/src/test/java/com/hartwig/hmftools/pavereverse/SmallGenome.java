package com.hartwig.hmftools.pavereverse;

import java.nio.file.Path;

import com.hartwig.hmftools.common.genome.tiny.ChromosomeSnippet;
import com.hartwig.hmftools.common.genome.tiny.TinyGenome;

public class SmallGenome extends TinyGenome
{
    public SmallGenome()
    {
        super();
        add(new ChromosomeSnippet("chr1", 10_000_000, 3_000_000, snippetPath( "chr1_part.txt")));
        add(new ChromosomeSnippet("chr1", 26_000_000, 1_000_000, snippetPath( "chr1_part_26.txt")));
        add(new ChromosomeSnippet("chr3", 10_000_000, 3_000_000, snippetPath( "chr3_part_10.txt")));
        add(new ChromosomeSnippet("chr4", 54_000_000, 1_000_000, snippetPath( "chr4_part_54.txt")));
        add(new ChromosomeSnippet("chr4", 105_000_000, 3_000_000, snippetPath( "chr4_part_105.txt")));
        add(new ChromosomeSnippet("chr5", 1_000_000, 1_000_000, snippetPath( "chr5_part_1.txt")));
        add(new ChromosomeSnippet("chr5", 68_000_000, 1_000_000, snippetPath( "chr5_part_68.txt")));
        add(new ChromosomeSnippet("chr7", 4_000_000, 1_000_000, snippetPath( "chr7_part_4.txt")));
        add(new ChromosomeSnippet("chr7", 55_000_000, 3_000_000, snippetPath( "chr7_part_55.txt")));
        add(new ChromosomeSnippet("chr7", 140_000_000, 10_000_000, snippetPath( "chr7_part140.txt")));
        add(new ChromosomeSnippet("chr17", 43_000_000, 1_000_000, snippetPath( "chr17_part_43.txt")));
        add(new ChromosomeSnippet("chr21", 37_000_000, 1_000_000, snippetPath( "chr21_part_37.txt")));
    }

    private Path snippetPath(String snippet)
    {
        return Path.of("src", "test", "resources", "tinygenome", snippet);
    }
}

