package com.hartwig.hmftools.gripss

data class BreakEnd(val startOrientation: Byte)
    : BreakJunction()

data class BreakPoint(val startOrientation: Byte, val chromosome: String, val position: Int, val endOrientation: Byte)
    : BreakJunction()

sealed class BreakJunction() {

    companion object Factory {
        private const val BREAKPOINT_REGEX = "^(.*)([\\[\\]])(.+)[\\[\\]](.*)\$"

        fun create(alt: String): BreakJunction {
            if (alt.startsWith(".")) {
                return BreakEnd(-1)
            }

            if (alt.endsWith(".")) {
                return BreakEnd(1)
            }

            return breakpoint(alt);
        }

        private fun breakpoint(alt: String): BreakJunction {
            val match = Regex(BREAKPOINT_REGEX).find(alt)!!
            val (initialSequence, bracket, location, finalSequence) = match.destructured
            val (chromosome, position) = location.split(":")

            val startOrientation: Byte = if (initialSequence.isNotEmpty()) 1 else -1
            val endOrientation: Byte = if (bracket == "]") 1 else -1

            return BreakPoint(startOrientation, chromosome, position.toInt(), endOrientation);
        }
    }
}
