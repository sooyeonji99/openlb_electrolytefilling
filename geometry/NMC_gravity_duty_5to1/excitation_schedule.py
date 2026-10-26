"""Expected acceleration for validating saved on/rest results."""
import math
import sys


def active(config, time):
    on = float(config['excitation_on_s'])
    rest = float(config['excitation_rest_s'])
    period = on + rest
    tolerance = 64 * sys.float_info.epsilon * max(period, abs(time))
    if time < 0 or time >= float(config['excitation_end_s']) - tolerance:
        return False
    local = time % period
    if local < tolerance or period - local < tolerance:
        local = 0
    return rest == 0 or local < on - tolerance


def excitation_at(config, time):
    if not active(config, time):
        return 0.0
    return (float(config['gravity_scale']) * float(config['g_reference_m_s2'])
            * math.sin(phase_at(config, time)))


def active_time(config, time):
    time = max(0., min(time, float(config['excitation_end_s'])))
    on = float(config['excitation_on_s'])
    period = on + float(config['excitation_rest_s'])
    tolerance = 64 * sys.float_info.epsilon * max(period, abs(time))
    sets = math.floor((time + tolerance) / period)
    return sets * on + min(on, max(0., time - sets * period))


def frequency_at(config, time):
    low = float(config['excitation_frequency_min_hz'])
    high = float(config['excitation_frequency_max_hz'])
    period = float(config['excitation_sweep_period_s'])
    return low + (high-low) * (1-math.cos(2*math.pi*time/period)) / 2


def phase_at(config, time):
    low = float(config['excitation_frequency_min_hz'])
    high = float(config['excitation_frequency_max_hz'])
    period = float(config['excitation_sweep_period_s'])
    return 2*math.pi*(high+low)/2*time - (high-low)/2*period*math.sin(2*math.pi*time/period)


def active_cycles(config, time):
    """Integrate frequency over on windows, excluding rest."""
    end = max(0., min(time, float(config['excitation_end_s'])))
    on = float(config['excitation_on_s'])
    period = on + float(config['excitation_rest_s'])
    total = 0.
    for k in range(math.ceil(end / period)):
        start = k * period
        stop = min(start+on, end)
        total += (phase_at(config, stop)-phase_at(config, start))/(2*math.pi)
    return total
