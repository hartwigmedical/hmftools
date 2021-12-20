package com.hartwig.hmftools.serve.transvar;

import java.io.FileNotFoundException;

import com.google.common.io.Resources;
import com.hartwig.hmftools.common.genome.refgenome.RefGenomeVersion;
import com.hartwig.hmftools.serve.EnsemblDataCacheTestFactory;

import org.jetbrains.annotations.NotNull;

public final class TransvarTestFactory {

    private static final String REF_GENOME_FASTA_FILE_37 = Resources.getResource("refgenome/v37/ref.fasta").getPath();
    private static final String REF_GENOME_FASTA_FILE_38 = Resources.getResource("refgenome/v38/ref.fasta").getPath();

    private TransvarTestFactory() {
    }

    @NotNull
    static Transvar testTransvar37(@NotNull TransvarProcess process) {
        return new Transvar(process, testInterpreter37(), EnsemblDataCacheTestFactory.create37());
    }

    @NotNull
    static TransvarInterpreter testInterpreter37() {
        try {
            return TransvarInterpreter.withRefGenome(RefGenomeVersion.V37, REF_GENOME_FASTA_FILE_37);
        } catch (FileNotFoundException exception) {
            throw new IllegalStateException("Cannot create test interpreter! Message=" + exception.getMessage());
        }
    }

    @NotNull
    static TransvarInterpreter testInterpreter38() {
        try {
            return TransvarInterpreter.withRefGenome(RefGenomeVersion.V38, REF_GENOME_FASTA_FILE_38);
        } catch (FileNotFoundException exception) {
            throw new IllegalStateException("Cannot create test interpreter! Message=" + exception.getMessage());
        }
    }
}
