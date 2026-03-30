

/*
Copyright (c) 2026 Sam78887
Licensed under SNSL v1.0
Created by: Sam78887
Non-Commercial Use Only
*/



/**

Atomic Lab Architect JAVA - ALAJ
 * Atomic Architect
 * A Java-based chemical analysis tool and molecular dynamics simulator.
 * * Features
 * * Electron Configuration
 * Automatically generates subshell notation (e.g., 1s^2 2s^2) for elements 1 through 118.
 * * VSEPR Theory
 * Predicts molecular geometry, bond angles, and steric numbers to determine the 3D shape.
 * * MD Simulation
 * Runs a digital lab simulation using a Lennard-Jones potential to model particle interactions.
 * * Physical State Prediction
 * Estimates boiling points and predicts whether a compound will be a solid, liquid, or gas.
 */
import java.util.*;

public class AtomicArchitect3 {

    static class Atom {
        String name, symbol, config;
        int valenceElectrons, valency, atomicNum;
        double weight, en, radius;

        static final String[] NAMES = {"", "Hydrogen", "Helium", "Lithium", "Beryllium", "Boron", "Carbon", "Nitrogen", "Oxygen", "Fluorine", "Neon", "Sodium", "Magnesium", "Aluminum", "Silicon", "Phosphorus", "Sulfur", "Chlorine", "Argon", "Potassium", "Calcium", "Scandium", "Titanium", "Vanadium", "Chromium", "Manganese", "Iron", "Cobalt", "Nickel", "Copper", "Zinc", "Gallium", "Germanium", "Assenic", "Selenium", "Bromine", "Krypton", "Rubidium", "Strontium", "Yttrium", "Zirconium", "Niobium", "Molybdenum", "Technetium", "Ruthenium", "Rhodium", "Palladium", "Silver", "Cadmium", "Indium", "Tin", "Antimony", "Tellurium", "Iodine", "Xenon", "Cesium", "Barium", "Lanthanum", "Cerium", "Praseodymium", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};
        static final String[] SYMBOLS = {"", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"};
        static final double[] WEIGHTS = new double[119];
        static final double[] EN_SCALE = new double[119];
        static final double[] RADII = new double[119];

static {
    // 1. Initialize with defaults to prevent nulls
    Arrays.fill(WEIGHTS, 0.0);
    Arrays.fill(EN_SCALE, 0.0);
    Arrays.fill(RADII, 0.0);

    // 2. Full Atomic Weights (Standard Atomic Masses)
    double[] allWeights = {
        0, 1.008, 4.0026, 6.94, 9.0122, 10.81, 12.011, 14.007, 15.999, 18.998, 20.18, // 0-10
        22.99, 24.305, 26.982, 28.085, 30.974, 32.06, 35.45, 39.948, 39.098, 40.078, // 11-20
        44.956, 47.867, 50.942, 51.996, 54.938, 55.845, 58.933, 58.693, 63.546, 65.38, // 21-30
        69.723, 72.63, 74.922, 78.971, 79.904, 83.798, 85.468, 87.62, 88.906, 91.224, // 31-40
        92.906, 95.95, 98, 101.07, 102.91, 106.42, 107.87, 112.41, 114.82, 118.71, // 41-50
        121.76, 127.6, 126.9, 131.29, 132.91, 137.33, 138.91, 140.12, 140.91, 144.24, // 51-60
        145, 150.36, 151.96, 157.25, 158.93, 162.5, 164.93, 167.26, 168.93, 173.05, // 61-70
        174.97, 178.49, 180.95, 183.84, 186.21, 190.23, 192.22, 195.08, 196.97, 200.59, // 71-80
        204.38, 207.2, 208.98, 209, 210, 222, 223, 226, 227, 232.04, 231.04, 238.03, // 81-92
        237, 244, 243, 247, 247, 251, 252, 257, 258, 259, 262, 267, 268, 271, 270, // 93-108
        270, 278, 281, 281, 285, 286, 289, 289, 293, 294, 294 // 109-118
    };
    System.arraycopy(allWeights, 0, WEIGHTS, 0, allWeights.length);

    // 3. Pauling Electronegativity (EN_SCALE)
    // Common elements manually set, others follow a calculated periodic trend
    double[] en = {
        0, 2.20, 0, 0.98, 1.57, 2.04, 2.55, 3.04, 3.44, 3.98, 0, // 0-10
        0.93, 1.31, 1.61, 1.90, 2.19, 2.58, 3.16, 0, 0.82, 1.00 // 11-20
    };
    System.arraycopy(en, 0, EN_SCALE, 0, en.length);
    // Fill remaining EN based on general trends (Group 17 high, Group 1 low)
    for(int i=21; i<=118; i++) {
        if (EN_SCALE[i] == 0) {
            if (i >= 57 && i <= 71) EN_SCALE[i] = 1.1; // Lanthanides
            else if (i >= 89 && i <= 103) EN_SCALE[i] = 1.3; // Actinides
            else EN_SCALE[i] = 1.5 + (0.01 * (i % 18)); // Rough heuristic for heavy metals
        }
    }
    EN_SCALE[47]=1.93; EN_SCALE[79]=2.54; EN_SCALE[17]=3.16; EN_SCALE[35]=2.96; EN_SCALE[53]=2.66;

    // 4. Atomic Radii (pm)
    // Simplified Radii trends: Decreases across period, increases down group
    for (int i = 1; i <= 118; i++) {
        int shell = (i <= 2) ? 1 : (i <= 10) ? 2 : (i <= 18) ? 3 : (i <= 36) ? 4 : (i <= 54) ? 5 : (i <= 86) ? 6 : 7;
        RADII[i] = 30 + (shell * 25) - ((i % 10) * 2);
    }
    // Specific Overrides for accuracy
    RADII[1]=37; RADII[6]=77; RADII[7]=75; RADII[8]=73; RADII[9]=71; 
    RADII[17]=99; RADII[47]=144; RADII[79]=144; RADII[80]=151;
}

        Atom(int num) {
            this.atomicNum = (num >= 1 && num <= 118) ? num : 0;
            this.name = NAMES[this.atomicNum];
            this.symbol = SYMBOLS[this.atomicNum];
            this.weight = WEIGHTS[this.atomicNum];
            this.en = EN_SCALE[this.atomicNum];
            this.radius = RADII[this.atomicNum];
            generateConfigAndValence(this.atomicNum);
        }

private void generateConfigAndValence(int n) {
    String[] subshells = {"1s", "2s", "2p", "3s", "3p", "4s", "3d", "4p", "5s", "4d", "5p", "6s", "4f", "5d", "6p", "7s", "5f", "6d", "7p"};
    int[] caps = {2, 2, 6, 2, 6, 2, 10, 6, 2, 10, 6, 2, 14, 10, 6, 2, 14, 10, 6};
    StringBuilder sb = new StringBuilder();
    int temp = n;

    // 1. Standard Aufbau Filling Loop
    for (int i = 0; i < subshells.length && temp > 0; i++) {
        int fill = Math.min(temp, caps[i]);
        sb.append(subshells[i]).append("^").append(fill).append(" ");
        temp -= fill;
    }
    this.config = sb.toString().trim();

    // 2. Exception Handler (Stability of Half/Full d-shells)
    // These override the string generated by the loop above
    if (n == 24) this.config = "1s^2 2s^2 2p^6 3s^2 3p^6 4s^1 3d^5";      // Chromium (Cr)
    else if (n == 29) this.config = "1s^2 2s^2 2p^6 3s^2 3p^6 4s^1 3d^10";     // Copper (Cu)
    else if (n == 42) this.config = "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^5";  // Molybdenum (Mo)
    else if (n == 47) this.config = "1s^2 2s^2 2p^6 3s^2 3p^6 4s^2 3d^10 4p^6 5s^1 4d^10"; // Silver (Ag)
    else if (n == 79) this.config = "[Xe] 4f^14 5d^10 6s^1"; // Gold (Au)

    // 3. Valence and Valency Logic
    if (n == 2 || n == 10 || n == 18 || n == 36 || n == 54 || n == 86 || n == 118) {
        this.valency = 0; 
        this.valenceElectrons = (n == 2) ? 2 : 8; // Correction: Helium has 2 valence electrons
    } else {
        int mod = n % 18;
        this.valenceElectrons = (mod == 0) ? 8 : (mod > 10 ? mod - 10 : (mod <= 2 ? mod : (mod <= 12 ? 2 : mod - 10)));
        
        // Manual Overrides for specific group behaviors
        if (n == 1) valenceElectrons = 1; 
        if (n == 6) valenceElectrons = 4;
        if (n == 7) valenceElectrons = 5;
        if (n == 8) valenceElectrons = 6;
        if (n == 9) valenceElectrons = 7;
        
        // Correcting Valence Electrons for the transition metals with exceptions
        if (n == 24 || n == 29 || n == 47 || n == 79) valenceElectrons = 1;

        this.valency = (valenceElectrons <= 4) ? valenceElectrons : 8 - valenceElectrons;
    }
}

    }
    static class Particle {
        double x, y, z, vx, vy, vz, fx, fy, fz;
        Particle(double x, double y, double z) { this.x = x; this.y = y; this.z = z; }
    }
    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        String choice = "y";

        while (choice.equalsIgnoreCase("y")) {
            System.out.print("\nHow many unique element types? ");
            int typeCount = sc.nextInt();
            Map<String, Integer> counts = new LinkedHashMap<>();
            List<Atom> atomList = new ArrayList<>();
            double totalWeight = 0;
            int totalValencySum = 0;
            int maxVal = -1;
            Atom central = null, secondary = null;
            int totalAtoms = 0;

            for (int i = 0; i < typeCount; i++) {
                System.out.print("Enter Atomic Number (1-118): ");
                int z = sc.nextInt();
                System.out.print("Enter Quantity: ");
                int qty = sc.nextInt();
                Atom cur = new Atom(z);
                System.out.println("-> " + cur.name + " [" + cur.config + "]");
                totalWeight += (cur.weight * qty);
                totalValencySum += (cur.valency * qty);
                totalAtoms += qty;
                counts.put(cur.symbol, qty);
                for(int j=0; j<qty; j++) atomList.add(cur);

                if (cur.valency > maxVal || (cur.valency == maxVal && qty == 1)) {
                    maxVal = cur.valency;
                    secondary = central;
                    central = cur;
                } else if (secondary == null && !cur.symbol.equals(central.symbol)) {
                    secondary = cur;
                }
            }

            System.out.println("\n--- Discovery Scanner Analysis ---");
            boolean valCheck = (totalValencySum % 2 == 0);
            boolean saturated = true;
            if (counts.containsKey("C")) {
                int others = totalAtoms - counts.get("C");
                if (others < (counts.get("C") / 2 + 1) && !counts.containsKey("O")) saturated = false;
            }
            

// --- UPDATED DISCOVERY SCANNER LOGIC ---
            if (totalAtoms == 1) {
                // Correctly identifies Silver (Ag), Gold (Au), or Helium (He)
                System.out.println("Status: [PURE ELEMENTAL LATTICE / MONATOMIC]");
            } 
            else if (!valCheck || totalValencySum < (2 * maxVal)) {
                // Catches radicals (like OH) or chemically impossible structures (like CH5)
                System.out.println("Status: [IMPOSSIBLE / UNSTABLE RADICAL]");
            } 
            else if (!saturated) {
                // Catches unsaturated hydrocarbons that might be highly reactive
                System.out.println("Status: [HIGHLY REACTIVE / THEORETICAL]");
            } 
            else {
                // Checks if peripheral atoms are physically too large for the central atom
                if (secondary != null && secondary.radius > central.radius * 1.5 && counts.get(secondary.symbol) >= 4) {
                    System.out.println("Status: [PHYSICALLY STRAINED]");
                } else {
                    System.out.println("Status: [STABLE COMPOUND]");
                }
            }

            System.out.println("\n--- Mass Composition ---");
            for (String s : counts.keySet()) {
                double atomW = 0;
                for (int i=1; i<119; i++) if (Atom.SYMBOLS[i].equals(s)) atomW = Atom.WEIGHTS[i];
                System.out.printf("%s: %.2f%%\n", s, (atomW * counts.get(s) / totalWeight) * 100);
            }

            predictPhysicalState(counts, totalWeight);

            System.out.println("\n--- Structural Details ---");
            System.out.println("Molecular Formula: " + generateFormula(counts));
            System.out.println("SMILES Notation: " + generateSMILES(central, counts));
            System.out.printf("Total Molecular Weight: %.5f g/mol\n", totalWeight);

            // Inside the while loop in main
System.out.print("\nAdd a second molecule for a reaction test? (y/n): ");
if (sc.next().equalsIgnoreCase("y")) {
    Map<String, Integer> countsB = new HashMap<>();
    System.out.print("Second Molecule - Unique element types: ");
    int typesB = sc.nextInt();
    for (int j = 0; j < typesB; j++) {
        System.out.print("Atomic Num (1-118): ");
        int zB = sc.nextInt();
        System.out.print("Quantity: ");
        int qB = sc.nextInt();
        countsB.put(Atom.SYMBOLS[zB], qB);
    }
    predictReaction(counts, countsB);
}
            
            
            if (secondary != null) {
                double diff = Math.abs(central.en - secondary.en);
                String bType = (diff >= 1.7) ? "Ionic" : (diff >= 0.5) ? "Polar Covalent" : "Non-Polar Covalent";
                System.out.println("Primary Bond Type: " + bType + " (Delta EN: " + String.format("%.2f", diff) + ")");
                
                // New: Spectroscopy & Bond Analysis
                double bondLength = (central.radius + secondary.radius) - 9 * diff;
                System.out.printf("Predicted Bond Length (%s-%s): %.2f pm\n", central.symbol, secondary.symbol, bondLength);
                double mu = (central.weight * secondary.weight) / (central.weight + secondary.weight);

                // --- Enhanced Bond & Spectroscopy Analysis ---
if (counts.containsKey("O") && counts.containsKey("H")) {
    Atom oxy = new Atom(8);
    Atom hyd = new Atom(1);
    
    System.out.println("\n--- Hydroxyl (-OH) Bond Analysis ---");
    
    // 1. Bond Polarity
    double diffOH = Math.abs(oxy.en - hyd.en);
    System.out.printf("O-H Delta EN: %.2f (Highly Polar Covalent)\n", diffOH);

    // 2. Bond Length (Schomaker-Stevenson)
    double bondLengthOH = (oxy.radius + hyd.radius) - 9 * diffOH;
    System.out.printf("Predicted O-H Length: %.2f pm\n", bondLengthOH);

    // 3. Vibrational Frequency (The "IR Stretch")
    double muAmu = (oxy.weight * hyd.weight) / (oxy.weight + hyd.weight);
    double muKg = muAmu * 1.660539e-27;
    
    // Force constant for O-H is typically much higher (~780 N/m)
    double kOH = 780.0; 
    double freqHz = (1.0 / (2.0 * Math.PI)) * Math.sqrt(kOH / muKg);
    double freqTHz = freqHz / 1.0e12;

    System.out.printf("O-H Stretching Frequency: %.2f THz (IR Active)\n", freqTHz);
    
    // 4. Hydrogen Bonding Check
    if (counts.get("O") >= 1 && counts.get("H") >= 1) {
        System.out.println("Intermolecular Force: Strong Hydrogen Bonding detected.");
    }
}


                //double vibFreq = (1.0 / (2 * Math.PI)) * Math.sqrt((500.0 * (central.en * secondary.en / 4.0)) / (mu * 1.66e-27)) * 1e-12 / 1e8;
                //System.out.printf("Vibration Frequency: ~%.2f THz\n", vibFreq);
                if (secondary != null) {
    // 1. Calculate Reduced Mass (mu) in amu
    double muAmu = (central.weight * secondary.weight) / (central.weight + secondary.weight);
    
    // 2. Convert mu to kg (CRITICAL STEP)
    double muKg = muAmu * 1.660539e-27;
    
    // 3. Estimate Force Constant (k) in N/m 
    // We scale based on Electronegativity as a proxy for bond strength
    double k = 500.0 * (1.0 + Math.abs(central.en - secondary.en) / 2.0);
    
    // 4. Calculate Frequency in Hz, then convert to THz
    // 1 THz = 10^12 Hz
    double freqHz = (1.0 / (2.0 * Math.PI)) * Math.sqrt(k / muKg);
    double freqTHz = freqHz / 1.0e12;

    System.out.printf("Vibration Frequency: %.2f THz\n", freqTHz);
}

/** 
            if (totalAtoms > 1) {
                int attached = totalAtoms - counts.get(central.symbol); 
                int lp = (central.valenceElectrons - (attached)) / 2;
                if (lp < 0) lp = 0;
                int sn = attached + lp;
                System.out.println("Central Atom: " + central.name);
                System.out.println("Steric Number: " + sn + " (Atoms: " + attached + ", Lone Pairs: " + lp + ")");
                predictGeometry(sn, lp);
            }
*/
            if (totalAtoms >= 1) {
    System.out.println("\n--- Multi-Center Geometry Analysis ---");
    
    // 1. Carbon Center Analysis
    int cCount = counts.getOrDefault("C", 0);
    if (cCount > 0) {
        System.out.println("At Carbon Center: Tetrahedral (109.5°)");
    }

    // 2. Oxygen (Hydroxyl) Center Analysis
    if (counts.containsKey("O") && counts.containsKey("H")) {
        // Oxygen in -OH has 2 bonds (C-O and O-H) and 2 lone pairs
        int snO = 4; 
        int lpO = 2;
        System.out.print("At Oxygen (-OH) Center: ");
        predictGeometry(snO, lpO); 
    }
            }
            System.out.print("\nRun Digital Lab Simulation (MD) for this molecule? (y/n): ");
            if (sc.next().equalsIgnoreCase("y")) {
                runLabSimulation(counts, totalWeight);
            }
            
            System.out.print("\nAnalyze another? (y/n): ");
            choice = sc.next();
            
    


            }}}


    public static void runLabSimulation(Map<String, Integer> counts, double molWeight) {
        List<Particle> particles = new ArrayList<>();
        double sigma = 3.5, epsilon = 0.2, targetTemp = 200.0;
        String molName = "Generic Compound";

        // Inside runLabSimulation
int cCount = counts.getOrDefault("C", 0);
boolean hBond = counts.containsKey("H") && (counts.containsKey("O") || counts.containsKey("N"));

// Sync the Lab Thermostat with the Physical State Prediction
double bpCelsius = -170.0 + (cCount * 32.0) + (molWeight * 0.2) + (hBond ? 160.0 : 0.0);
if (cCount == 0 && counts.getOrDefault("O", 0) == 1) bpCelsius = 100.0; // Water override

targetTemp = bpCelsius + 273.15; // Convert to Kelvin for the MD loop

// Calibration: Adjust "Stickiness" (Epsilon) based on chain length
// This prevents the "Freezing" effect in the small 40-molecule box
epsilon = hBond ? 0.55 : (cCount > 4 ? 0.22 : 0.15);


        // Heuristic Intelligence Layer
        boolean hBonding = counts.containsKey("H") && (counts.containsKey("N") || counts.containsKey("O") || counts.containsKey("F"));
        if (counts.getOrDefault("C",0)==1 && counts.getOrDefault("H",0)==4) { epsilon=0.148; sigma=3.73; targetTemp=111.6; molName="Methane (CH4)"; }
        else if (counts.getOrDefault("O",0)==1 && counts.getOrDefault("H",0)==2) { epsilon=0.650; sigma=3.16; targetTemp=373.15; molName="Water (H2O)"; }
        else {

            
            // 1. Get structural counts for the formula
            
            int hCount = counts.getOrDefault("H", 0);
            int oCount = counts.getOrDefault("O", 0);
            int nCount = counts.getOrDefault("N", 0);
            int sCount = counts.getOrDefault("S", 0);

            // 2. Linear Surface Area Logic: +22.0K per Carbon fixes the "Heptanol Gap"
            // Base value (40.0) + Carbon chain influence + heavy atom mass influence
            double surfaceAreaBase = (cCount * 22.0) + (molWeight * 0.4) + 40.0;

            // 3. Polarity Logic: Differentiates OH (Alcohol) vs NH (Amine) vs SH (Thiol)
            double polarityBonus = 0;
            if (hCount > 0) {
                if (oCount > 0) polarityBonus = 165.0;      // Strong H-Bonding
                else if (nCount > 0) polarityBonus = 85.0;  // Medium H-Bonding
                else if (sCount > 0) polarityBonus = 20.0;  // Very Weak (Sulfur)
            }

            // 4. Final Target Calculation
            targetTemp = surfaceAreaBase + polarityBonus + 200.0;
            
            // Adjust Epsilon (Stickiness) based on polarity
            epsilon = (oCount > 0 || nCount > 0) ? 0.55 : 0.25; 
            
            molName = generateFormula(counts);
        }
            /** 
            //epsilon = hBonding ? 0.45 : 0.18;
            //targetTemp = (molWeight * 1.2) + (hBonding ? 180 : 60);
            // Enhanced Universal Heuristic
double baseBP = 120.0 * Math.log10(molWeight + 10.0); // Mass-based scaling
double hBondBonus = hBonding ? 150.0 : 0.0;           // Polarity jump
targetTemp = baseBP + hBondBonus;

// Ensure we don't drop below a reasonable 'liquid' floor for organics
if (hBonding && targetTemp < 330) targetTemp = 330; 

molName = generateFormula(counts); // Dynamically name it
        } */

        System.out.println("\n--- STABILIZED DIGITAL LAB: " + molName + " ---");
        int numMols = 40; double dt = 0.5;
        int side = (int) Math.ceil(Math.pow(numMols, 1.0/3.0));
        for (int i = 0; i < numMols; i++) particles.add(new Particle((i%side)*6.0, ((i/side)%side)*6.0, (i/(side*side))*6.0));

        for (int step = 0; step < 1001; step++) {
            for (Particle p : particles) { p.fx = p.fy = p.fz = 0; }
            for (int i = 0; i < particles.size(); i++) {
                for (int j = i + 1; j < particles.size(); j++) {
                    Particle p1 = particles.get(i), p2 = particles.get(j);
                    double dx = p1.x-p2.x, dy = p1.y-p2.y, dz = p1.z-p2.z;
                    double r2 = dx*dx + dy*dy + dz*dz + 0.5;
                    double r = Math.sqrt(r2);
                    if (r < 12.0) {
                        double r6 = Math.pow(sigma / r, 6);
                        double fS = (48 * epsilon / r2) * (r6 * r6 - 0.5 * r6);
                        if (fS > 80) fS = 80;
                        p1.fx += fS * dx; p1.fy += fS * dy; p1.fz += fS * dz;
                        p2.fx -= fS * dx; p2.fy -= fS * dy; p2.fz -= fS * dz;
                    }
                }
            }
            double currentKE = 0;
            for (Particle p : particles) {
                double acc = 0.0004184 / molWeight;
                p.vx += p.fx * acc * dt; p.vy += p.fy * acc * dt; p.vz += p.fz * acc * dt;
                p.vx = Math.max(-1.5, Math.min(1.5, p.vx));
                p.x += p.vx * dt; p.y += p.vy * dt; p.z += p.vz * dt;
                currentKE += (p.vx*p.vx + p.vy*p.vy + p.vz*p.vz);
            }
            double currentT = (currentKE * molWeight) / (3 * numMols * 0.008314);
            if (step % 10 == 0) {
                double scale = Math.sqrt(1.0 + 0.1 * (targetTemp / (currentT + 0.1) - 1.0));
                for (Particle p : particles) { p.vx *= scale; p.vy *= scale; p.vz *= scale; }
            }
            if (step % 200 == 0) System.out.printf("Step %d | Temp: %.2f K | Status: %s\n", step, currentT, (currentT >= targetTemp - 5 ? "STABILIZED" : "HEATING"));
        }
    }

public static void predictPhysicalState(Map<String, Integer> counts, double weight) {
    int cCount = counts.getOrDefault("C", 0);
    int hCount = counts.getOrDefault("H", 0);
    int oCount = counts.getOrDefault("O", 0);
    int nCount = counts.getOrDefault("N", 0);
    
    // Identify Substance Type
    boolean isMetal = false;
    boolean isNobleGas = false;
    String mainSymbol = counts.size() == 1 ? counts.keySet().iterator().next() : "";

    if (counts.size() == 1) {
        // Noble Gas Check
        if ("He".equals(mainSymbol) || "Ne".equals(mainSymbol) || "Ar".equals(mainSymbol) || 
            "Kr".equals(mainSymbol) || "Xe".equals(mainSymbol) || "Rn".equals(mainSymbol)) {
            isNobleGas = true;
        }
        // Metallic Check (Simplified: Elements with weight > 40 that aren't Non-metals)
        else if (weight > 40 && !mainSymbol.equals("I") && !mainSymbol.equals("Br") && !mainSymbol.equals("S")) {
            isMetal = true;
        }
    }

    double bpCelsius;

    if (isMetal) {
        // Metallic bonding is extremely strong (Sea of Electrons)
        bpCelsius = 1500.0 + (weight * 5.0); 
    } else if (isNobleGas) {
        // Extremely weak Van der Waals forces
        bpCelsius = -272.0 + (weight * 0.4);
    } else if (cCount == 0 && oCount == 1 && hCount == 2) {
        bpCelsius = 100.0; // Water override
    } else {
        // Original Organic/Covalent logic
        boolean hBond = hCount > 0 && (oCount > 0 || nCount > 0);
        double baseBP = -170.0; 
        double chainEffect = (cCount * 32.0); 
        double massEffect = (weight * 0.22);
        double polarityEffect = hBond ? 160.0 : 0.0;
        if (cCount > 0 && oCount > 0 && hCount > 0) polarityEffect += 40.0;

        bpCelsius = baseBP + chainEffect + massEffect + polarityEffect;
    }

    System.out.println("\n--- Physical State Prediction ---");
    System.out.printf("Estimated Boiling Point: ~%.1f°C\n", bpCelsius);
    
    if (isMetal) {
        System.out.println("Predicted Phase: Solid (Metallic Lattice)");
    } else if (bpCelsius < 25) {
        System.out.println("Predicted Phase: Gas (at 25°C)");
    } else if (cCount > 12 || weight > 300) {
        System.out.println("Predicted Phase: Solid (Molecular/Waxy)");
    } else {
        System.out.println("Predicted Phase: Liquid (at 25°C)");
    }
}

public static String generateFormula(Map<String, Integer> counts) {
    // Specialized format for Alcohols (e.g., C2H5OH instead of C2H6O)
    if (counts.getOrDefault("C", 0) > 0 && counts.getOrDefault("O", 0) > 0 && counts.getOrDefault("H", 0) > 1) {
        int c = counts.get("C");
        int h = counts.get("H");
        return String.format("C%sH%sOH", (c > 1 ? c : ""), (h - 1 > 1 ? (h - 1) : ""));
    }

    // Hill System: Carbon first, then Hydrogen, then Alphabetical
    StringBuilder sb = new StringBuilder();
    List<String> sorted = new ArrayList<>(counts.keySet());
    Collections.sort(sorted, (a, b) -> {
        if (a.equals("C")) return -1;
        if (b.equals("C")) return 1;
        if (a.equals("H")) return -1;
        if (b.equals("H")) return 1;
        return a.compareTo(b);
    });

    for (String s : sorted) {
        int qty = counts.get(s);
        sb.append(s).append(qty > 1 ? qty : "");
    }
    return sb.toString();
}



public static String generateSMILES(Atom central, Map<String, Integer> counts) {
    StringBuilder sm = new StringBuilder();

    // 1. Handle Organic Chains (C-C-C...)
    if (counts.containsKey("C")) {
        for (int i = 0; i < counts.get("C"); i++) sm.append("C");
    } else {
        // For Inorganics, start with the Central Atom
        sm.append(central.symbol);
    }

    // 2. Add Functional Groups (O, N, S, P, etc.)
    // We skip Hydrogen in SMILES because it's usually "implicit"
    for (Map.Entry<String, Integer> entry : counts.entrySet()) {
        String sym = entry.getKey();
        if (!sym.equals("C") && !sym.equals("H") && !sym.equals(central.symbol)) {
            sm.append("(").append(sym).append(")");
        }
    }
    
    // 3. Special case for Alcohols: If Oxygen is present, ensure it shows as the -OH group
    if (counts.containsKey("O") && !sm.toString().contains("O")) {
        sm.append("O");
    }

    return sm.toString();
}
/** 
    public static void predictGeometry(int sn, int lp) {
        System.out.print("Molecular Geometry: ");
        if (sn == 2) System.out.println("Linear (180 degrees)");
        else if (sn == 3) System.out.println(lp == 0 ? "Trigonal Planar (120 degrees)" : "Bent (~118 degrees)");
        else if (sn == 4) {
            if (lp == 0) System.out.println("Tetrahedral (109.5 degrees)");
            else if (lp == 1) System.out.println("Trigonal Pyramidal (~107 degrees)");
            else System.out.println("Bent (~104.5 degrees)");
        } else System.out.println("Complex/Advanced Geometry");
    }
*/

    // Updated predictGeometry method:
public static void predictGeometry(int sn, int lp) {
    if (sn == 2) System.out.println("Linear (180°)");
    else if (sn == 3) System.out.println(lp == 0 ? "Trigonal Planar (120°)" : "Bent (~118°)");
    else if (sn == 4) {
        if (lp == 0) System.out.println("Tetrahedral (109.5°)");
        else if (lp == 1) System.out.println("Trigonal Pyramidal (~107°)");
        else if (lp == 2) System.out.println("Bent (~104.5°) due to Lone Pair Repulsion");
    } else {
        System.out.println("Complex/Advanced Geometry");
    }
}
/* 
public static void predictReaction(Map<String, Integer> molA, Map<String, Integer> molB) {
    System.out.println("\n--- ⚗️ Reaction Prediction Lab ---");
    
    int cA = molA.getOrDefault("C", 0), hA = molA.getOrDefault("H", 0), oA = molA.getOrDefault("O", 0);
    int cB = molB.getOrDefault("C", 0), hB = molB.getOrDefault("H", 0), oB = molB.getOrDefault("O", 0);

    // 1. Combustion Check (Fuel + Oxygen)
    if (cA > 0 && oB >= 2 && cB == 0) {
        System.out.println("Reaction: [COMBUSTION / OXIDATION]");
        double enthalpy = (cA * 393.5) + (hA * 141.8); // kJ/mol
        System.out.printf("Products: %dCO2 + %dH2O\n", cA, hA/2);
        System.out.printf("Energy Released: -%.2f kJ/mol (Exothermic)\n", enthalpy);
    } 
    // 2. Esterification (Alcohol + Carboxylic Acid)
    else if (oA >= 1 && hA >= 1 && cB >= 1 && oB >= 2) {
        System.out.println("Reaction: [FISCHER ESTERIFICATION]");
        System.out.println("Mechanism: Nucleophilic Acyl Substitution");
        



        System.out.printf("Product: Ester (C%dH%dO2) + Water\n", (cA + cB), (hA + hB - 2));
        System.out.println("Note: Likely produces a fruity or floral aroma.");
    } 
    else {
        System.out.println("Status: No spontaneous reaction detected at STP.");
    }
        
}
*/

public static void predictReaction(Map<String, Integer> molA, Map<String, Integer> molB) {
    System.out.println("\n--- ⚗️ Reaction Prediction Lab ---");
    
    // Molecule A Profile
    int cA = molA.getOrDefault("C", 0), hA = molA.getOrDefault("H", 0), oA = molA.getOrDefault("O", 0);
    int nA = molA.getOrDefault("N", 0), sA = molA.getOrDefault("S", 0), fA = molA.getOrDefault("F", 0);
    
    // Molecule B Profile
    int cB = molB.getOrDefault("C", 0), hB = molB.getOrDefault("H", 0), oB = molB.getOrDefault("O", 0);
    int clB = molB.getOrDefault("Cl", 0), sB = molB.getOrDefault("S", 0), fB = molB.getOrDefault("F", 0);

    // 1. Combustion Check (Fuel + Oxygen)
    if (cA > 0 && oB >= 2 && cB == 0) {
        System.out.println("Reaction: [COMBUSTION / OXIDATION]");
        double enthalpy = (cA * 393.5) + (hA * 141.8); 
        System.out.printf("Products: %dCO2 + %dH2O\n", cA, hA/2);
        System.out.printf("Energy Released: -%.2f kJ/mol (Highly Exothermic)\n", enthalpy);
    } 
    

    // 2. Fischer Esterification (Alcohol + Carboxylic Acid)
    else if(oA >= 1 && hA >= 1 && cB >= 1 && oB >= 2 && hB >= 1) {
        System.out.println("Reaction: [FISCHER ESTERIFICATION]");
        System.out.println("Mechanism: Nucleophilic Acyl Substitution");
        int totalC = cA + cB;
        int totalH = hA + hB - 2;
        System.out.printf("Product: Ester (C%dH%dO2) + Water (H2O)\n", totalC, totalH);
        
        // Aroma Database Logic
        if (totalC == 7) System.out.println("Aroma Detected: Banana (Isoamyl Acetate)");
        else if (totalC == 9) System.out.println("Aroma Detected: Orange (Octyl Acetate)");
        else System.out.println("Note: Likely produces a fruity or floral aroma.");
    }

    // 3. Acid-Base Neutralization (Amine + Acid)
    else if (nA >= 1 && clB >= 1) {
        System.out.println("Reaction: [ACID-BASE NEUTRALIZATION]");
        System.out.println("Mechanism: Proton Transfer (Brønsted–Lowry)");
        System.out.println("Product: Ammonium Chloride (NH4Cl) - White Salt/Smoke.");
    }

    // 4. Halogenation (Organic + Fluorine/Chlorine)
    else if (cA > 0 && (fB >= 2 || clB >= 2)) {
        System.out.println("Reaction: [HALOGENATION]");
        System.out.println("Mechanism: Free Radical Substitution");
        System.out.println("Status: Carbon-Halogen bonds forming. Potential Teflon precursor.");
    }

    // 5. Synthesis of Sulfides (Metal/H2 + Sulfur)
    else if (sB > 0 && (hA > 0 || molA.containsKey("Fe") || molA.containsKey("Cu"))) {
        System.out.println("Reaction: [SULFIDE FORMATION]");
        if (hA > 0) System.out.println("Product: Hydrogen Sulfide (H2S) - Warning: Toxic/Rotten Egg Odor.");
        else System.out.println("Product: Metallic Sulfide (Tarnish/Mineral Layer).");
    }

    // 6. Precipitation (Silver + Halogen) - Specific Test
    else if (molA.containsKey("Ag") && clB > 0) {
        System.out.println("Reaction: [PRECIPITATION]");
        System.out.println("Product: Silver Chloride (AgCl) - White insoluble solid.");
    }

    else {
        System.out.println("Status: No spontaneous reaction detected at STP.");
        System.out.println("Suggestion: Try increasing Digital Lab Temperature to overcome Activation Energy.");
    }
}
}

