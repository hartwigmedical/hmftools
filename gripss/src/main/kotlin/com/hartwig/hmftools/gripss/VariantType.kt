package com.hartwig.hmftools.gripss

import kotlin.math.abs

data class Single(
        override val insertSequence: String,
        val startOrientation: Byte)
    : VariantType()

data class Translocation(
        override val insertSequence: String,
        override val otherChromosome: String,
        override val otherPosition: Int,
        override val startOrientation: Byte,
        override val endOrientation: Byte)
    : Paired()

data class Inversion(
        override val insertSequence: String,
        override val otherChromosome: String,
        override val otherPosition: Int,
        override val startOrientation: Byte,
        override val endOrientation: Byte)
    : Paired()

data class Deletion(
        override val insertSequence: String,
        override val otherChromosome: String,
        override val otherPosition: Int,
        override val startOrientation: Byte,
        override val endOrientation: Byte,
        val length: Int)
    : Paired()

data class Duplication(
        override val insertSequence: String,
        override val otherChromosome: String,
        override val otherPosition: Int,
        override val startOrientation: Byte,
        override val endOrientation: Byte,
        val length: Int) : Paired()

data class Insertion(
        override val insertSequence: String,
        override val otherChromosome: String,
        override val otherPosition: Int,
        override val startOrientation: Byte,
        override val endOrientation: Byte,
        val length: Int) : Paired()

sealed class Paired : VariantType() {
    abstract val otherChromosome: String
    abstract val otherPosition: Int
    abstract val startOrientation: Byte
    abstract val endOrientation: Byte
}

sealed class VariantType {
    abstract val insertSequence: String

    companion object Factory {
        private const val BREAKPOINT_REGEX = "^(.*)([\\[\\]])(.+)[\\[\\]](.*)\$"

        fun create(chromosome: String, position: Int, ref: String, alt: String): VariantType {
            if (alt.startsWith(".")) {
                return Single(alt.substring(ref.length, alt.length - 1), -1)
            }

            if (alt.endsWith(".")) {
                return Single(alt.substring(1, alt.length - ref.length), 1)
            }

            return paired(chromosome, position, ref, alt);
        }

        private fun paired(chromosome: String, position: Int, ref: String, alt: String): Paired {
            val match = Regex(BREAKPOINT_REGEX).find(alt)!!
            val (initialSequence, bracket, location, finalSequence) = match.destructured
            val locationList = location.split(":");
            val (otherChromosome, otherPosition) = Pair(locationList[0], locationList[1].toInt())

            val endOrientation: Byte = if (bracket == "]") 1 else -1

            val startOrientation: Byte
            val insertSequence: String
            if (initialSequence.isNotEmpty()) {
                startOrientation = 1
                insertSequence = initialSequence.substring(ref.length, initialSequence.length)
            } else {
                startOrientation = -1
                insertSequence = finalSequence.substring(0, finalSequence.length - ref.length)
            }

            if (chromosome != otherChromosome) {
                return Translocation(insertSequence, otherChromosome, otherPosition, startOrientation, endOrientation)
            }

            if (startOrientation == endOrientation) {
                return Inversion(insertSequence, otherChromosome, otherPosition, startOrientation, endOrientation)
            }

            val length = abs(position - otherPosition)

            if ((position <= otherPosition && startOrientation == (-1).toByte()) || (position >= otherPosition && startOrientation == (1).toByte())) {
                return Duplication(insertSequence, otherChromosome, otherPosition, startOrientation, endOrientation, length)
            }

            return if (length <= 1)
                Insertion(insertSequence, otherChromosome, otherPosition, startOrientation, endOrientation, length)
            else
                Deletion(insertSequence, otherChromosome, otherPosition, startOrientation, endOrientation, length)
        }
    }
}
