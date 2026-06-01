from rational_design.ai_core import AIExpertAgent


class DummyBackend:
    pass


def test_ai_expert_prompt_keeps_full_pcr_advisory_role():
    agent = AIExpertAgent(DummyBackend())

    assert "KHÔNG chỉ là một công cụ thiết kế primer" in agent.system_instruction_vi
    assert "tối ưu wet-lab" in agent.system_instruction_vi
    assert "propose_validation" in agent.system_instruction_vi
    assert "propose_multiplex" in agent.system_instruction_vi
    assert "degenerate_primers" in agent.system_instruction_vi
    assert "IUPAC" in agent.system_instruction_vi

    assert "NOT only a primer-design executor" in agent.system_instruction_en
    assert "wet-lab optimization" in agent.system_instruction_en
    assert "propose_validation" in agent.system_instruction_en
    assert "propose_multiplex" in agent.system_instruction_en
    assert "degenerate_primers" in agent.system_instruction_en
    assert "IUPAC" in agent.system_instruction_en


def test_ai_expert_prompt_separates_advice_from_runnable_jobs():
    agent = AIExpertAgent(DummyBackend())

    assert "KHÔNG xuất block JSON" in agent.system_instruction_vi
    assert "DO NOT emit a runnable JSON block" in agent.system_instruction_en
