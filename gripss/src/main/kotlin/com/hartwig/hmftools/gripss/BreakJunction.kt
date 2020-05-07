package com.hartwig.hmftools.gripss

data class BreakEnd(override val insertSequence: String, val startOrientation: Byte)
    : BreakJunction()

data class BreakPoint(override val insertSequence: String, val startOrientation: Byte, val chromosome: String, val position: Int, val endOrientation: Byte)
    : BreakJunction()

sealed class BreakJunction() {
    abstract val insertSequence: String

    companion object Factory {
        private const val BREAKPOINT_REGEX = "^(.*)([\\[\\]])(.+)[\\[\\]](.*)\$"

        fun create(ref: String, alt: String): BreakJunction {
            if (alt.startsWith(".")) {
                return BreakEnd(alt.substring(ref.length, alt.length - 1), -1)
            }

            if (alt.endsWith(".")) {
                return BreakEnd(alt.substring(1, alt.length - ref.length), 1)
            }

            return breakpoint(ref, alt);
        }

        private fun breakpoint(ref: String, alt: String): BreakJunction {
            val match = Regex(BREAKPOINT_REGEX).find(alt)!!
            val (initialSequence, bracket, location, finalSequence) = match.destructured
            val (chromosome, position) = location.split(":")

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

            return BreakPoint(insertSequence, startOrientation, chromosome, position.toInt(), endOrientation);
        }
    }
}
