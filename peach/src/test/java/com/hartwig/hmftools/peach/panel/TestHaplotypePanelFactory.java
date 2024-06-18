package com.hartwig.hmftools.peach.panel;

import java.util.Map;

import com.google.common.collect.ImmutableList;
import com.hartwig.hmftools.common.genome.chromosome.HumanChromosome;
import com.hartwig.hmftools.peach.event.VariantHaplotypeEvent;
import com.hartwig.hmftools.peach.haplotype.DefaultHaplotype;
import com.hartwig.hmftools.peach.haplotype.NonDefaultHaplotype;

import org.jetbrains.annotations.NotNull;

public class TestHaplotypePanelFactory
{
    @NotNull
    public static HaplotypePanel createDefaultTestHaplotypePanel()
    {
        GeneHaplotypePanel dpydGeneHaplotypePanel = createDpydGeneHaplotypePanel();
        GeneHaplotypePanel ugt1a1GeneHaplotypePanel = createUgt1a1GeneHaplotypePanel();
        Map<String, GeneHaplotypePanel> geneToHaplotypePanel = Map.of(
                "DPYD", dpydGeneHaplotypePanel,
                "UGT1A1", ugt1a1GeneHaplotypePanel
        );
        return new HaplotypePanel(geneToHaplotypePanel);
    }

    @NotNull
    public static HaplotypePanel createDpydTestHaplotypePanel()
    {
        GeneHaplotypePanel dpydGeneHaplotypePanel = createDpydGeneHaplotypePanel();
        Map<String, GeneHaplotypePanel> geneToHaplotypePanel = Map.of(
                "DPYD", dpydGeneHaplotypePanel
        );
        return new HaplotypePanel(geneToHaplotypePanel);
    }

    @NotNull
    public static HaplotypePanel createUgt1a1TestHaplotypePanel()
    {
        GeneHaplotypePanel ugt1a1GeneHaplotypePanel = createUgt1a1GeneHaplotypePanel();
        Map<String, GeneHaplotypePanel> geneToHaplotypePanel = Map.of(
                "UGT1A1", ugt1a1GeneHaplotypePanel
        );
        return new HaplotypePanel(geneToHaplotypePanel);
    }

    @NotNull
    private static GeneHaplotypePanel createDpydGeneHaplotypePanel()
    {
        DefaultHaplotype dpydDefaultHaplotype = new DefaultHaplotype("*9A", false, ImmutableList.of());
        ImmutableList<NonDefaultHaplotype> dpydNonDefaultHaplotypes = ImmutableList.of(
                new NonDefaultHaplotype("*1", true,
                        ImmutableList.of(new VariantHaplotypeEvent(HumanChromosome._1, 98348885, "G", "A"))
                ),
                new NonDefaultHaplotype("*2A", false,
                        ImmutableList.of(
                                new VariantHaplotypeEvent(HumanChromosome._1, 97915614, "C", "T"),
                                new VariantHaplotypeEvent(HumanChromosome._1, 98348885, "G", "A")
                        )
                ),
                new NonDefaultHaplotype("*7", false,
                        ImmutableList.of(
                                new VariantHaplotypeEvent(HumanChromosome._1, 98205966, "GATGA", "G"),
                                new VariantHaplotypeEvent(HumanChromosome._1, 98348885, "G", "A")
                        )
                ),
                new NonDefaultHaplotype("*13", false,
                        ImmutableList.of(
                                new VariantHaplotypeEvent(HumanChromosome._1, 97981343, "A", "C"),
                                new VariantHaplotypeEvent(HumanChromosome._1, 98348885, "G", "A")
                        )
                ),
                new NonDefaultHaplotype("*B3", false,
                        ImmutableList.of(
                                new VariantHaplotypeEvent(HumanChromosome._1, 98039419, "C", "T"),
                                new VariantHaplotypeEvent(HumanChromosome._1, 98045449, "G", "C"),
                                new VariantHaplotypeEvent(HumanChromosome._1, 98348885, "G", "A")
                        )
                ),
                new NonDefaultHaplotype("2846A>T", false,
                        ImmutableList.of(
                                new VariantHaplotypeEvent(HumanChromosome._1, 97547947, "T", "A"),
                                new VariantHaplotypeEvent(HumanChromosome._1, 98348885, "G", "A")
                        )
                )
        );
        return new GeneHaplotypePanel(dpydDefaultHaplotype, dpydNonDefaultHaplotypes, "*1");
    }

    @NotNull
    private static GeneHaplotypePanel createUgt1a1GeneHaplotypePanel()
    {
        DefaultHaplotype ugt1a1DefaultHaplotype = new DefaultHaplotype("*1", true, ImmutableList.of());
        ImmutableList<NonDefaultHaplotype> ugt1a1NonDefaultHaplotypes = ImmutableList.of(
                new NonDefaultHaplotype("*28", false,
                        ImmutableList.of(new VariantHaplotypeEvent(HumanChromosome._2, 234668879, "C", "CAT"))
                ),
                new NonDefaultHaplotype("*37", false,
                        ImmutableList.of(new VariantHaplotypeEvent(HumanChromosome._2, 234668879, "C", "CATAT"))
                ),
                new NonDefaultHaplotype("*36", false,
                        ImmutableList.of(new VariantHaplotypeEvent(HumanChromosome._2, 234668879, "CAT", "C"))
                ),
                new NonDefaultHaplotype("*6", false,
                        ImmutableList.of(new VariantHaplotypeEvent(HumanChromosome._2, 234669144, "G", "A"))
                )
        );
        return new GeneHaplotypePanel(ugt1a1DefaultHaplotype, ugt1a1NonDefaultHaplotypes, "*1");
    }
}
