import {
  BarChart,
  Callout,
  Card,
  CardBody,
  CardHeader,
  CollapsibleSection,
  Grid,
  H1,
  H2,
  H3,
  Pill,
  Row,
  Stack,
  Stat,
  Table,
  Text,
  TodoListCard,
  UsageBar,
  useHostTheme,
} from "cursor/canvas";

const M_TARGET = 41.284;
const I_STAR = 1.2993;

const terminalRoutes = [
  { name: "Phase D Route D", omega: 80.72, iTerm: 112.84, status: "refuted" as const },
  { name: "Phase E Route D", omega: 68.69, iTerm: 78.8, status: "refuted" as const },
  { name: "L-0073 equal_12", omega: 54.23, iTerm: 37.92, status: "best" as const },
  { name: "fine_low", omega: 55.01, iTerm: 40.13, status: "worse" as const },
  { name: "equal_24", omega: 63.4, iTerm: 63.85, status: "worse" as const },
  { name: "equal_48 (Route D)", omega: 65.75, iTerm: 70.5, status: "worse" as const },
];

const lowSlabRoutes = [
  { name: "triad_min_defect (all-IC sketch)", omega: 51.0, allIc: true },
  { name: "L-0073 terminal", omega: 54.23, allIc: true },
  { name: "weighted_min_defect", omega: 57.3, allIc: true },
  { name: "sos_greedy one-pol", omega: 21.71, allIc: false },
  { name: "sos_band123", omega: 17.54, allIc: false },
  { name: "M target", omega: M_TARGET, allIc: true },
];

const lemmas = [
  { id: "L-0072", status: "route_refuted", note: "Per-shell Phase E" },
  { id: "L-0073", status: "best_known_bound", note: "Terminal equal_12" },
  { id: "L-0074", status: "partition_refuted", note: "Ladder 24/48/fine_low" },
  { id: "L-0075", status: "low_slab_hybrid", note: "Sketch all-IC Omega~51" },
];

const solutions = [
  {
    id: "A1",
    title: "L-0027 all-IC C <= C_dagger",
    effort: "Math + empirical",
    impact: "High",
    closes: "All-IC if proved",
  },
  {
    id: "A2",
    title: "High-slab L-0024/L-0026 rigor",
    effort: "Medium",
    impact: "High",
    closes: "Partial",
  },
  {
    id: "A3",
    title: "New terminal structure (non-partition)",
    effort: "High research",
    impact: "Unknown",
    closes: "Possible",
  },
  {
    id: "B1",
    title: "Weighted triad refinement",
    effort: "Hours-days",
    impact: "Medium",
    closes: "Open",
  },
  {
    id: "B2",
    title: "SOS band extension + rigor",
    effort: "Days",
    impact: "Subclass",
    closes: "Subclass only",
  },
  {
    id: "C1",
    title: "More clusters / Route D",
    effort: "—",
    impact: "None",
    closes: "Refuted L-0074",
  },
];

const artifacts = [
  "certificates/CERT-L0073-C0008-route-a-prime-best.json",
  "experiments/terminal_weighted/shell_manifest_cluster_best_full.json",
  "experiments/terminal_weighted/route_a_ladder_full.json",
  "experiments/terminal_weighted/hybrid_low_slab_summary.json",
  "reports/C0008_MASTER_HANDOFF.md",
  "run_route_a_ladder_resilient.cmd",
  "reproduce_c0008_certificate.ps1",
];

const todoNext = [
  { id: "1", content: "Prove all-IC C <= C_dagger (L-0027)", status: "pending" as const },
  { id: "2", content: "Rigor high-slab + low-slab combination", status: "pending" as const },
  { id: "3", content: "Do NOT promote C-0008 to proved", status: "completed" as const },
  { id: "4", content: "Terminal partition family exhausted (L-0074)", status: "completed" as const },
];

export default function C0008MasterHandoff() {
  const { tokens } = useHostTheme();
  const chartData = terminalRoutes.map((r) => ({
    label: r.name.replace("L-0073 ", ""),
    value: r.omega,
  }));

  return (
    <Stack gap={24} style={{ padding: 24, maxWidth: 960 }}>
      <Stack gap={8}>
        <H1>C-0008 Master Handoff</H1>
        <Text tone="secondary">
          Galerkin T³ N≤24 · full dealias · status exploring · target M≈41.284
        </Text>
        <Row gap={12}>
          <Stat label="Best terminal Ω" value="54.23" tone="warning" />
          <Stat label="Target M" value="41.28" tone="default" />
          <Stat label="Gap Ω" value="+12.9" tone="error" />
          <Stat label="Best I_hi" value="37.92" tone="warning" />
          <Stat label="I_* lo" value="1.30" tone="default" />
        </Row>
      </Stack>

      <Callout tone="warning">
        Verifier FAIL is expected. Do not mark C-0008 proved. Terminal cluster
        partition family refuted (L-0074). Best bound: L-0073 equal_12.
      </Callout>

      <Card>
        <CardHeader title="Terminal routes — Ω(T)^hi vs M" />
        <CardBody>
          <BarChart
            data={chartData}
            series={[{ dataKey: "value", name: "Ω(T)^hi", tone: "accent" }]}
            referenceLines={[{ value: M_TARGET, label: "M target", tone: "error" }]}
            xLabel="Route"
            yLabel="Ω(T)^hi"
            caption="Source: Phase D/E, L-0073, Route A ladder full · Aug 2026"
          />
        </CardBody>
      </Card>

      <Grid columns={2} gap={16}>
        <Card>
          <CardHeader title="Low-slab hybrid (L-0075)" />
          <CardBody>
            <Table
              columns={[
                { key: "name", header: "Route" },
                { key: "omega", header: "Ω worst", align: "right" },
                {
                  key: "allIc",
                  header: "All-IC",
                  render: (v: boolean) => (
                    <Pill tone={v ? "success" : "warning"} size="sm">
                      {v ? "sketch" : "subclass"}
                    </Pill>
                  ),
                },
              ]}
              data={lowSlabRoutes}
            />
            <Text size="sm" tone="secondary" style={{ marginTop: 12 }}>
              Best all-IC sketch: triad_min_defect Ω≈51.0 (gap +9.7 vs M). Subclass
              SOS Ω≈17.5 does not close all-IC.
            </Text>
          </CardBody>
        </Card>

        <Card>
          <CardHeader title="Lemmas registry" />
          <CardBody>
            <Stack gap={8}>
              {lemmas.map((l) => (
                <Row key={l.id} gap={8} style={{ alignItems: "center" }}>
                  <Pill tone="accent" size="sm">
                    {l.id}
                  </Pill>
                  <Text weight="medium">{l.status}</Text>
                  <Text tone="secondary" size="sm">
                    {l.note}
                  </Text>
                </Row>
              ))}
            </Stack>
          </CardBody>
        </Card>
      </Grid>

      <CollapsibleSection
        header={
          <Row gap={8}>
            <H3>Compute timeline</H3>
            <Pill size="sm">~30h+ total</Pill>
          </Row>
        }
      >
        <Stack gap={8}>
          <UsageBar
            segments={[
              { label: "Phase E", value: 24, color: "category1" },
              { label: "Route A ladder", value: 5, color: "category2" },
              { label: "L-0073", value: 1.2, color: "category3" },
              { label: "Audits", value: 0.1, color: "category4" },
            ]}
            caption="Wall-clock hours (approximate)"
          />
          <Text size="sm" tone="secondary">
            Phase D ~40min · L-0073 ~70min · Ladder 4 configs ~5h · Phase E ~24h+
          </Text>
        </Stack>
      </CollapsibleSection>

      <Card>
        <CardHeader title="Possible solutions (prioritized)" />
        <CardBody>
          <Table
            columns={[
              { key: "id", header: "ID", width: 48 },
              { key: "title", header: "Route" },
              { key: "effort", header: "Effort" },
              { key: "impact", header: "Impact" },
              { key: "closes", header: "Closes?" },
            ]}
            data={solutions}
          />
        </CardBody>
      </Card>

      <Grid columns={2} gap={16}>
        <Card>
          <CardHeader title="Checkpoint / resume" />
          <CardBody>
            <Stack gap={6}>
              <Text size="sm">run_route_a_ladder_resilient.cmd — per config ckpt</Text>
              <Text size="sm">run_upgrade_best_resilient.ps1 — per shell ckpt</Text>
              <Text size="sm">Cluster D² pass: NO intra-chunk ckpt (~70 min/config)</Text>
              <Text size="sm">SOS cluster D=1772: ~113h, no Python ckpt</Text>
            </Stack>
          </CardBody>
        </Card>

        <Card>
          <CardHeader title="Key artifacts" />
          <CardBody>
            <Stack gap={4}>
              {artifacts.map((a) => (
                <Text key={a} size="sm" style={{ fontFamily: tokens.font.mono }}>
                  {a}
                </Text>
              ))}
            </Stack>
          </CardBody>
        </Card>
      </Grid>

      <TodoListCard todos={todoNext} defaultExpanded />

      <Text size="sm" tone="secondary">
        Full detail: reports/C0008_MASTER_HANDOFF.md · Reproduce: reproduce_c0008_certificate.ps1
      </Text>
    </Stack>
  );
}
