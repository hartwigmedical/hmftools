
plugins {
    id("com.hartwig.java-conventions")
}

description = "HMF Tools - TELO"

dependencies {
    implementation(project(":hmf-common"))
    implementation(libs.tablesaw.core)
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
    testImplementation(libs.junit)
}

shadowJar {
    mainClass.set("com.hartwig.hmftools.telo.TeloApplication")
}
