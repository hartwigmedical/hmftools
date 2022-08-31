package com.hartwig.hmftools.cider

import org.eclipse.collections.api.collection.ImmutableCollection
import org.eclipse.collections.api.set.SetIterable

// we use immutable collections here, data can be accessed by multiple threads
interface CiderGeneDatastore
{
    fun getAnchorSequenceSet(geneType: VJGeneType): SetIterable<String>
    fun getByAnchorSequence(anchorSeq: String): ImmutableCollection<VJAnchorTemplate>
    fun getByAnchorSequence(geneType: VJGeneType, anchorSeq: String): ImmutableCollection<VJAnchorTemplate>
    fun getByAnchorGeneLocation(vjAnchorReferenceLocation: VJAnchorReferenceLocation): ImmutableCollection<VJAnchorTemplate>
    fun getVJAnchorReferenceLocations(): SetIterable<VJAnchorReferenceLocation>
    fun getIgConstantRegions(): ImmutableCollection<IgTcrConstantRegion>
}